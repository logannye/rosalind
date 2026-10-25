"""Small, bounded contracts shared by the producer and saved-only verifier."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import threading
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
UINT64_MAX = 2**64 - 1
COUNTS = ('callable_depth', 'alt_count', 'alt_fwd', 'alt_rev')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    if Path(path).stat().st_size > 32 << 20:
        raise ValueError('JSON metadata exceeds 32 MiB')
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Atomic replacement is only used inside a newly owned run directory.
    temporary = path.with_name(path.name + '.partial')
    with temporary.open('x') as output:
        json.dump(value, output, sort_keys=True, indent=2, allow_nan=False)
        output.write('\n')
    temporary.replace(path)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def safe_path(root, name):
    if not isinstance(name, str) or '\\' in name or not name:
        raise ValueError('invalid bundle-relative path')
    relative = PurePosixPath(name)
    if relative.is_absolute() or any(p in ('', '.', '..') for p in name.split('/')):
        raise ValueError('invalid bundle-relative path')
    root = Path(root).resolve()
    candidate = root
    for part in relative.parts:
        candidate = candidate / part
        if candidate.is_symlink():
            raise ValueError('bundle symlink refused')
    if not candidate.resolve().is_relative_to(root):
        raise ValueError('path escapes bundle')
    return candidate


def inventory(root, exclude=()):
    root = Path(root)
    files = {}
    for path in sorted(root.rglob('*')):
        if path.is_symlink():
            raise ValueError('symlink in retained directory')
        if path.is_file():
            relative = path.relative_to(root).as_posix()
            if relative not in exclude:
                files[relative] = {'sha256': sha256(path), 'bytes': path.stat().st_size}
    return files


def verify_inventory(root, entries):
    if not isinstance(entries, dict) or not entries or len(entries) > 10000:
        raise ValueError('invalid bounded inventory')
    for name, record in entries.items():
        path = safe_path(root, name)
        if not path.is_file() or path.stat().st_size != record['bytes'] or sha256(path) != record['sha256']:
            raise ValueError(f'inventory mismatch: {name}')


def size(root):
    return sum(record['bytes'] for record in inventory(root).values())


def load_protocol(path):
    if sha256(path) != sha256(HERE / 'protocol-v1.json'):
        raise ValueError('protocol differs from this exact source revision')
    return read_json(path)


def probes(path):
    rows = []
    with Path(path).open() as source:
        for line in source:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 8:
                raise ValueError('expected eight-column probe VCF')
            name, pos, _, ref, alt, *_ = fields
            if len(ref) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or len(alt) != 1 or ref == alt:
                raise ValueError('invalid probe SNV')
            rows.append(dict(contig=name, pos=int(pos), ref=ref, alt=alt))
    return rows


def write_sites(path, rows, contigs):
    with Path(path).open('x') as output:
        output.write('##fileformat=VCFv4.2\n')
        for item in contigs:
            output.write(f"##contig=<ID={item['name']},length={item['length']}>\n")
        output.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for row in rows:
            output.write(f"{row['contig']}\t{row['pos']}\t.\t{row['ref']}\t{row['alt']}\t.\t.\t.\n")


def request_definitions(rows, contigs):
    ranks = {row['name']: rank for rank, row in enumerate(contigs)}
    lengths = {row['name']: row['length'] for row in contigs}
    keys = [(ranks[r['contig']], r['pos']) for r in rows]
    if len(rows) != 1000 or keys != sorted(set(keys)):
        raise ValueError('exactly 1,000 unique dictionary-ordered probes required')
    if any(r['contig'] != 'chr20' or not 10000001 <= r['pos'] <= 11000000
           or r['pos'] > lengths[r['contig']] or r['alt'] != next(b for b in 'ACGT' if b != r['ref']) for r in rows):
        raise ValueError('probe workload or original ALT differs from frozen definition')
    q2 = [dict(r, alt=[b for b in 'ACGT' if b != r['ref']][1]) for r in rows[::2]]
    return [dict(id=name, sites=f'requests/{name}.vcf', rows=chosen,
                 thresholds=thresholds, expected=f'expected/{name}.json')
            for name, chosen, thresholds in [('Q1', rows[:250], []), ('Q2', q2, []), ('Q3', rows, [5, 10, 20])]]


def freeze_queries(probes_vcf, destination, contigs):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    rows = probes(probes_vcf)
    requests = request_definitions(rows, contigs)
    (destination / 'probes.vcf').write_bytes(Path(probes_vcf).read_bytes())
    for request in requests:
        write_sites(destination / (request['id'] + '.vcf'), request['rows'], contigs)
    # Choose missing coordinate solely from the mask, before reading any counts.
    stored = {r['pos'] for r in rows}
    absent = next(p for p in range(10000001, 11000001) if p not in stored)
    (destination / 'absent.bed').write_text(f'chr20\t{absent-1}\t{absent}\n')
    wrong = dict(rows[0])
    wrong['ref'] = next(b for b in 'ACGT' if b not in (wrong['ref'], wrong['alt']))
    write_sites(destination / 'wrong-ref.vcf', [wrong], contigs)
    write_json(destination / 'freeze.json', {'requests': requests, 'files': inventory(destination)})
    return requests


def native_table(path):
    with Path(path).open() as source:
        reader = csv.DictReader(source, delimiter='\t')
        if not reader.fieldnames or reader.fieldnames[0] != '#contig':
            raise ValueError('unexpected native TSV schema')
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError('malformed TSV row')
            row['contig'] = row.pop('#contig')
            yield row


def uint64(value):
    if isinstance(value, bool) or not str(value).isdigit():
        raise ValueError('expected uint64')
    parsed = int(value)
    if not 0 <= parsed <= UINT64_MAX:
        raise ValueError('uint64 overflow')
    return parsed


def reduce_native(path, request):
    """Reducer under test: strict one row per request, including observed zero."""
    expected = {(r['contig'], r['pos']): r for r in request['rows']}
    seen, result = set(), []
    for row in native_table(path):
        key = row['contig'], int(row['pos'])
        if key in seen or key not in expected:
            raise ValueError('unexpected or duplicate native locus')
        seen.add(key)
        r = expected[key]
        if row['ref'] != r['ref'] or row['requested_alts'] != r['alt']:
            raise ValueError('native REF/ALT mismatch')
        alt = r['alt'].lower()
        answer = dict(r, callable_depth=uint64(row['callable_depth']), alt_count=uint64(row[alt]),
                      alt_fwd=uint64(row[alt+'_fwd']), alt_rev=uint64(row[alt+'_rev']))
        for threshold in request['thresholds']:
            answer[f'depth_ge_{threshold}'] = answer['callable_depth'] >= threshold
        result.append(answer)
    if seen != set(expected) or [(r['contig'], r['pos']) for r in result] != list(expected):
        raise ValueError('missing or out-of-order native locus')
    return result


def oracle_answers(path, request):
    """Independent expectation path; does not call the native reducer/normalizer."""
    with Path(path).open() as source:
        reader = csv.reader(source, delimiter='\t')
        header = next(reader)
        column = {name: i for i, name in enumerate(header)}
        observations = {}
        for values in reader:
            key = values[column['#contig']], int(values[column['pos']])
            if key in observations:
                raise ValueError('duplicate oracle locus')
            observations[key] = values
    answers = []
    for selected in request['rows']:
        values = observations[(selected['contig'], selected['pos'])]
        if values[column['ref']] != selected['ref']:
            raise ValueError('oracle REF mismatch')
        base = selected['alt'].lower()
        depth = int(values[column['callable_depth']])
        answer = {'contig': selected['contig'], 'pos': selected['pos'], 'ref': selected['ref'],
                  'alt': selected['alt'], 'callable_depth': depth,
                  'alt_count': int(values[column[base]]), 'alt_fwd': int(values[column[base+'_fwd']]),
                  'alt_rev': int(values[column[base+'_rev']])}
        answer.update({f'depth_ge_{t}': depth >= t for t in request['thresholds']})
        if any(not 0 <= answer[field] <= UINT64_MAX for field in COUNTS):
            raise ValueError('invalid oracle integer')
        answers.append(answer)
    return answers


def command(argv, destination, expect_failure=False):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    argv = list(map(str, argv))
    started = time.perf_counter()
    with (destination / 'stdout.log').open('xb') as stdout, (destination / 'stderr.log').open('xb') as stderr:
        code = subprocess.run(argv, stdout=stdout, stderr=stderr, check=False).returncode
    record = {'argv': argv, 'exit_code': code, 'wall_seconds': time.perf_counter()-started,
              'stdout_sha256': sha256(destination/'stdout.log'), 'stderr_sha256': sha256(destination/'stderr.log')}
    write_json(destination/'command.json', record)
    if (code == 0) == expect_failure:
        raise RuntimeError(f'unexpected command exit {code}: {argv}; retained logs: {destination}')
    return record


def native_query(binary, manifest, request_path, output):
    return [str(binary), 'dataset', 'extract', '--dataset', str(manifest), '--sites', str(request_path),
            '--fields', 'depths,alleles,strands', '--format', 'tsv', '--memory-budget-mb', '256',
            '--output', str(output), '--manifest', str(output)+'.manifest.json']


def measure(argv, destination):
    """Cohort benchmark's 10 ms simultaneous RSS sampler plus parent-clock ACKs.

    Worker stdout is a control channel; native output is captured separately.
    Each ACK follows verified publication, never just native process exit.
    """
    import psutil
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    argv = list(map(str, argv))
    stopped = threading.Event()
    peak, samples, prefixes = None, 0, []
    started = time.perf_counter()
    with (destination/'stderr.log').open('xb') as stderr, (destination/'stdout.log').open('x') as log:
        process = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=stderr, text=True, bufsize=1)
        def sample():
            nonlocal peak, samples
            try:
                parent = psutil.Process(process.pid)
            except psutil.Error:
                return
            while not stopped.is_set():
                try:
                    rss, observed = 0, False
                    for item in [parent, *parent.children(recursive=True)]:
                        try:
                            rss += item.memory_info().rss
                            observed = True
                        except psutil.Error:
                            pass
                    if observed and rss > 0:
                        peak = max(peak or 0, rss)
                        samples += 1
                except psutil.Error:
                    pass
                stopped.wait(0.01)
        monitor = threading.Thread(target=sample, daemon=True)
        monitor.start()
        try:
            for line in process.stdout:
                observed_at = time.perf_counter()-started
                log.write(line)
                log.flush()
                try:
                    signal = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if 'completed_query' in signal:
                    if signal['completed_query'] != f'Q{len(prefixes)+1}':
                        process.kill()
                        raise ValueError('invalid worker completion sequence')
                    prefixes.append(dict(signal, parent_prefix_seconds=observed_at))
            code = process.wait()
            elapsed = time.perf_counter()-started
        finally:
            stopped.set()
            monitor.join()
            if process.poll() is None:
                process.kill()
                process.wait()
            process.stdout.close()
    record = {'argv': argv, 'exit_code': code, 'outer_wall_seconds': elapsed,
              'sampled_process_tree_peak_rss_bytes': peak, 'rss_samples': samples,
              'rss_sample_interval_ms': 10, 'prefixes': prefixes,
              'rss_scope': 'simultaneous worker + descendants; excludes monitoring driver; may miss peaks or double-count shared pages'}
    write_json(destination/'measurement.json', record)
    return record
