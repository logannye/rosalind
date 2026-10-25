#!/usr/bin/env python3
"""One-command fixed real-read case study. Never builds or installs anything."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import statistics
import sys
import time
import traceback

from case_common import (HERE, ROOT, command, freeze_queries, inventory, load_protocol, measure,
                         module, native_query, native_table, oracle_answers, read_json, reduce_native,
                         sha256, verify_inventory, write_json)


def logical_bytes(path):
    return sum(p.stat().st_size for p in Path(path).rglob('*') if p.is_file())


def representatives():
    # Existing source preparation, general-CIGAR oracle argv, and exact TSV comparison.
    sys.path.insert(0, str(ROOT / 'benchmarks/evidence'))
    # representative.py's historical `from run import ...` must resolve to its
    # sibling, even when this entry point is imported by a test as `run`.
    previous = sys.modules.get('run')
    try:
        module('run', ROOT / 'benchmarks/evidence/run.py')
        return module('case_representative', ROOT / 'benchmarks/evidence/representative.py')
    finally:
        if previous is None:
            sys.modules.pop('run', None)
        else:
            sys.modules['run'] = previous


def fresh_command(binary, state, selection, output, fields='depths,alleles,strands', cache=None, fmt='tsv'):
    rep = representatives()
    workload = dict(state['workload'])
    workload['selection'] = dict(workload['selection'], path=str(selection))
    argv = rep.native_command(binary, workload, state['protocol']['native'], state['protocol']['oracle'],
                              output, Path(str(output)+'.manifest.json'), cache=cache)
    argv[argv.index('--fields')+1] = fields
    argv[argv.index('--format')+1] = fmt
    return argv


def receipt_verify(binary, output, log):
    return command([binary, 'verify', '--manifest', str(output)+'.manifest.json', '--json'], log)


def arrow_oracle_check(arrow_path, oracle_path):
    import pyarrow as pa
    from streaming_pysam import HEADER
    expected = iter(native_table(oracle_path))
    compared = 0
    with pa.ipc.open_stream(arrow_path) as stream:
        for batch in stream:
            if batch.num_rows > 1024:
                raise ValueError('native Arrow batch exceeds 1024')
            for row in batch.to_pylist():
                normalized = {}
                for key, value in row.items():
                    if key.endswith('_histogram'):
                        normalized[key] = ','.join(f'{i}:{n}' for i, n in enumerate(value) if n) or '.'
                    elif key == 'requested_alts':
                        normalized[key] = ','.join(value) or '.'
                    else:
                        normalized[key] = str(value)
                other = next(expected, None)
                if normalized != other:
                    raise ValueError(f'full Arrow/oracle mismatch at row {compared}: {normalized} != {other}')
                compared += 1
    if next(expected, None) is not None or compared != 1000:
        raise ValueError('full Arrow/oracle cardinality mismatch')
    return {'loci': compared, 'fields': len(HEADER), 'mismatches': 0, 'all_legacy_fields': True}


def worker(args):
    state = read_json(args.state)
    root = Path(state['root'])
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    binary = args.binary
    if args.worker == 'oracle':
        record = command(representatives().oracle_command(state['workload'], state['protocol']['oracle'],
                         output/'stats.json'), output/'command')
        shutil.copy2(output/'command/stdout.log', root/'oracle.tsv')
        write_json(output/'result.json', record)
        return
    if args.worker == 'materialize':
        arrow = output/'evidence.arrow'
        records = [command(fresh_command(binary, state, root/'requests/probes.vcf', arrow,
                    fields='all', cache=output/'cache', fmt='arrow-ipc'), output/'extract')]
        records.append(receipt_verify(binary, arrow, output/'verify-output'))
        receipt = read_json(str(arrow)+'.manifest.json')
        leaf = Path(receipt['measurements']['execution.evidence_dataset_manifest']).parent
        shutil.copytree(leaf, root/'native')
        manifest = root/'native/evidence-dataset.manifest.json'
        records.append(command([binary, 'dataset', 'verify', '--dataset', manifest], output/'verify-dataset'))
        verified = read_json(output/'verify-dataset/stdout.log')
        if verified['verified_loci'] != 1000 or verified['original_sources_rehashed'] is not False:
            raise ValueError('unexpected saved dataset verification')
        full = output/'full.tsv'
        argv = native_query(binary, manifest, root/'requests/probes.vcf', full)
        argv[argv.index('--fields')+1] = 'all'
        records.append(command(argv, output/'full-tsv'))
        records.append(receipt_verify(binary, full, output/'verify-tsv'))
        comparison = representatives().compare(root/'oracle.tsv', full)
        write_json(output/'comparison.json', comparison)
        if not comparison['equal']:
            raise ValueError('native full legacy fields differ from independent raw oracle')
        arrow_comparison = arrow_oracle_check(arrow, root/'oracle.tsv')
        values = receipt['measurements']
        # Missing values remain unknown and fail required admission checks.
        observed = {key: values.get(key) for key in ('execution.worker_count', 'execution.microtile_bases',
                    'execution.microtiles', 'execution.record_visits', 'execution.plan')}
        if int(observed['execution.worker_count']) != 1 or not 1 <= int(observed['execution.microtile_bases']) <= 16384:
            raise ValueError('native admitted configuration differs from fixed protocol')
        write_json(output/'result.json', {'commands': records, 'full_tsv': comparison,
                    'full_arrow': arrow_comparison, 'observed': observed, 'receipt_measurements': values})
        return
    if args.worker == 'convert':
        import checked_parquet
        descriptor = read_json(root/'native/dataset.descriptor.json')
        semantics = {key: descriptor[key] for key in ('semantics', 'counting_unit', 'sampling', 'has_reference',
                     'profile', 'sample_scope', 'fields', 'schema_version', 'selection', 'sources')}
        semantics['protocol_sha256'] = state['protocol_sha256']
        converted = checked_parquet.convert_arrow(root/'materialize/evidence.arrow', root/'parquet',
                      semantics=semantics, contig_order=[r['name'] for r in state['contigs']], expected_rows=1000)
        write_json(output/'result.json', converted)
        return
    observations = []
    try:
        for request in state['requests']:
            current = output/request['id']
            current.mkdir()
            started = time.perf_counter()
            verification_seconds = 0.0
            metrics = {}
            if args.worker in ('A', 'B'):
                native = current/'evidence.tsv'
                if args.worker == 'A':
                    argv = fresh_command(binary, state, root/request['sites'], native)
                else:
                    argv = native_query(binary, root/state['native'], root/request['sites'], native)
                execution = command(argv, current/'extract')
                verification = receipt_verify(binary, native, current/'verify')
                verification_seconds += verification['wall_seconds']
                answer = reduce_native(native, request)
                receipt = read_json(str(native)+'.manifest.json')
                metrics = {'native_subprocesses': 2, 'receipt_measurements': receipt.get('measurements', {}),
                           'native_output_sha256': sha256(native), 'native_exit_code': execution['exit_code']}
                if args.worker == 'B':
                    values = metrics['receipt_measurements']
                    # Absent counters are unknown and cannot satisfy this gate.
                    if values.get('execution.alignment_record_visits') != '0' or values.get('original_sources_rehashed') != 'false':
                        raise ValueError('saved query source-absence counters missing or nonzero')
            else:
                import checked_parquet
                answer, metrics = checked_parquet.query(root/state['parquet'], request['rows'],
                    thresholds=request['thresholds'], expected_inventory_sha256=state['parquet_inventory_sha256'])
                metrics['native_subprocesses'] = 0
                verification_seconds += metrics.get('verification_seconds', 0.0)
            write_json(current/'answer.json', answer)
            check_started = time.perf_counter()
            expected = read_json(root/request['expected'])
            if answer != expected or read_json(current/'answer.json') != expected:
                raise ValueError(f'{args.worker}/{request["id"]} scientific mismatch')
            verification_seconds += time.perf_counter()-check_started
            record = {'path': args.worker, 'query': request['id'], 'status': 'completed',
                      'query_wall_seconds': time.perf_counter()-started, 'verification_seconds': verification_seconds,
                      'rows': len(answer), 'answer_sha256': sha256(current/'answer.json'), 'metrics': metrics}
            observations.append(record)
            write_json(output/'observations.json', observations)
            print(json.dumps({'completed_query': request['id'], 'answer_sha256': record['answer_sha256']}), flush=True)
    except Exception as error:
        observations.append({'path': args.worker, 'query': request['id'], 'status': 'failed', 'error': str(error)})
        for pending in state['requests'][len(observations):]:
            observations.append({'path': args.worker, 'query': pending['id'], 'status': 'unstarted'})
        write_json(output/'observations.json', observations)
        raise


def synthetic_control(binary, root):
    """Reuse authored specimen C, whose stored position 40 has no observations."""
    destination = root/'synthetic'
    command([sys.executable, ROOT/'examples/cohort-reanalysis/prepare.py', destination/'fixture'],
            root/'logs/synthetic-prepare')
    raw = destination/'fixture/raw'
    output = destination/'evidence.tsv'
    command([binary, 'analyze', 'evidence', '--reference', raw/'reference.fa', '--alignments', raw/'specimen-c.bam',
             '--regions', raw/'specimen-c.bed', '--sample', 'sample-c', '--fields', 'depths,alleles,strands',
             '--memory-budget-mb', '256', '--workers', '1', '--format', 'tsv', '--cache-dir', destination/'cache',
             '--output', output], root/'logs/synthetic-extract')
    receipt_verify(binary, output, root/'logs/synthetic-verify')
    rows = list(native_table(output))
    zero = next(row for row in rows if int(row['pos']) == 40)
    if int(zero['callable_depth']) != 0:
        raise ValueError('authored stored zero control failed')
    leaf = Path(read_json(str(output)+'.manifest.json')['measurements']['execution.evidence_dataset_manifest']).parent
    shutil.copytree(leaf, root/'synthetic-zero/native')
    (root/'synthetic-zero/zero.bed').write_text('synthetic\t39\t40\n')
    (root/'synthetic-zero/absent.bed').write_text('synthetic\t40\t41\n')
    shutil.copy2(ROOT/'examples/cohort-reanalysis/expected.json', root/'synthetic-zero/authored-expected.json')
    return {'label': 'authored synthetic specimen C, separate from HG002', 'position': 40, 'callable_depth': 0}


def summary(report):
    setup = report['stages']['materialize']['outer_wall_seconds']
    conversion = report['stages']['convert']['outer_wall_seconds']
    results = {}
    for path in ('A', 'B', 'D'):
        trials = [r for r in report['trials'] if r['path'] == path]
        query_stats, prefix_stats = {}, {}
        added = 0 if path == 'A' else setup + (conversion if path == 'D' else 0)
        for index, query in enumerate(('Q1', 'Q2', 'Q3')):
            times = [r['observations'][index]['query_wall_seconds'] for r in trials]
            checks = [r['observations'][index]['verification_seconds'] for r in trials]
            prefixes = [r['prefixes'][index]['parent_prefix_seconds'] for r in trials]
            query_stats[query] = {'median': statistics.median(times), 'min': min(times), 'max': max(times),
                                  'verification_median': statistics.median(checks)}
            prefix_stats[query] = {'observed_parent_seconds': prefixes, 'median': statistics.median(prefixes),
                'min': min(prefixes), 'max': max(prefixes), 'setup_added_seconds': added,
                'setup_inclusive_median': added+statistics.median(prefixes),
                'setup_inclusive_min': added+min(prefixes), 'setup_inclusive_max': added+max(prefixes)}
        peaks = [r['sampled_process_tree_peak_rss_bytes'] for r in trials]
        if path != 'A':
            peaks.append(report['stages']['materialize']['sampled_process_tree_peak_rss_bytes'])
        if path == 'D':
            peaks.append(report['stages']['convert']['sampled_process_tree_peak_rss_bytes'])
        results[path] = {'queries_seconds': query_stats, 'prefixes_seconds': prefix_stats,
                         'workflow_max_sampled_process_tree_rss_bytes': max(peaks) if all(p is not None for p in peaks) else None,
                         'whole_worker_seconds': [r['outer_wall_seconds'] for r in trials]}
    return results


def human_report(report):
    text = ['# HG002 reproducibility case study', '', f"Producer status: **{report['run_status']}**. Clean replay is a separate acceptance gate.",
            '', 'D shares native extraction with B; it is a downstream representation comparison.',
            'Setup is measured once. The totals below add that observation to each parent-clock prefix.', '',
            '| Path | Common setup (s) | Conversion (s) | Q1 / Q2 / Q3 median (s) | K=3 including setup (s) | Minimum cache bytes | Sampled workflow RSS bytes | Exactness |',
            '|---|---:|---:|---|---:|---:|---:|---|']
    for path, result in report.get('summary', {}).items():
        q = result['queries_seconds']
        setup = report['stages']['materialize']['outer_wall_seconds'] if path != 'A' else 0
        conversion = report['stages']['convert']['outer_wall_seconds'] if path == 'D' else 0
        retained = report['storage_logical_bytes'].get({'A':'no_evidence_cache','B':'native_minimum_replay','D':'parquet_minimum_replay'}[path], 0)
        text.append(f"| {path} | {setup:.6f} | {conversion:.6f} | {q['Q1']['median']:.6f} / {q['Q2']['median']:.6f} / {q['Q3']['median']:.6f} | {result['prefixes_seconds']['Q3']['setup_inclusive_median']:.6f} | {retained} | {result['workflow_max_sampled_process_tree_rss_bytes']} | 9/9 oracle-equal |")
    text += ['', f"Acquisition plus preparation: {report['stages']['acquisition_preparation']['outer_wall_seconds']:.6f} s; {report['acquired_bytes']:,} newly acquired bytes, {report['authenticated_source_bytes']:,} authenticated source bytes.",
             f"Independent raw oracle: {report['stages']['oracle']['outer_wall_seconds']:.6f} s (evaluation cost, outside path totals).",
             'Common setup includes extraction, receipt verification, portable copy, full dataset verification, all-field TSV export/verification and exact TSV/Arrow-oracle checks.',
             'D conversion includes full-schema readback and inventory creation. A has no evidence cache but requires the retained raw BAM/reference inputs.',
             'Minimum cache bytes exclude the complete diagnostic bundle and query output receipts.', '',
             '| Path / prefix | Observed parent time median [range], s | Setup-inclusive median [range], s |',
             '|---|---|---|']
    for path, result in report.get('summary', {}).items():
        for query, p in result['prefixes_seconds'].items():
            text.append(f"| {path} / {query} | {p['median']:.6f} [{p['min']:.6f}, {p['max']:.6f}] | {p['setup_inclusive_median']:.6f} [{p['setup_inclusive_min']:.6f}, {p['setup_inclusive_max']:.6f}] |")
    for path in ('B','D'):
        crosses = [q for q in ('Q1','Q2','Q3') if report['summary'][path]['prefixes_seconds'][q]['setup_inclusive_median'] < report['summary']['A']['prefixes_seconds'][q]['setup_inclusive_median']]
        text.append(f"{path} setup-inclusive median below A at: {', '.join(crosses) if crosses else 'no observed prefix through three queries'}. This describes these measured prefixes only.")
    text += ['', 'All individual times, ranges, prefixes, exits, hashes and resource samples are in report.json.',
             'Acquisition, oracle evaluation, controls, source rechecks and bundle publication are separate stages.',
             'No cold-cache, p99, hard memory limit, clinical value, or outside-user adoption is established.',
             'Separate verification times are already in the totals. Native integrated hashing is not separately isolated; D explicitly times inventory and payload integrity checks.',
             'The earlier cohort analyzer no-go and frozen HG001 practical-benefit failure remain in force.', '']
    return '\n'.join(text)


def produce(args):
    import psutil
    from build_identity import authenticate
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    report = {'schema_version': 1, 'run_status': 'running', 'stages': {}, 'trials': [], 'failures': [],
              'clean_replay_status': 'unstarted', 'stopping_rule': 'one narrowly justified correction maximum; no timing selection'}
    report['stage_status'] = {name: 'unstarted' for name in ('identity', 'small_preflight', 'acquisition_preparation',
        'query_freeze', 'oracle', 'materialize', 'convert', 'trials', 'controls', 'source_rechecks', 'bundle')}
    try:
        protocol = load_protocol(args.protocol)
        identity = authenticate(args.binary, args.build_identity, args.expected_source)
        report['stage_status']['identity'] = 'completed'
        report['identity'] = identity
        report['protocol_sha256'] = sha256(args.protocol)
        shutil.copy2(args.protocol, root/'protocol-v1.json')
        shutil.copy2(args.build_identity, root/'build-identity.json')
        free = shutil.disk_usage(root).free
        available = psutil.virtual_memory().available
        report['resource_preflight'] = {'free_scratch_bytes': free, 'available_ram_bytes': available}
        if free < protocol['minimum_free_scratch_bytes'] or available < protocol['minimum_available_ram_bytes']:
            raise ValueError('insufficient resources before public-data acquisition')
        def stage(name, argv):
            report['stage_status'][name] = 'running'
            write_json(root/'report.json', report)
            result = measure(argv, root/'logs'/name)
            report['stages'][name] = result
            write_json(root/'report.json', report)
            if result['exit_code'] != 0:
                raise RuntimeError(f'{name} failed; retained logs and partial files')
            if result['sampled_process_tree_peak_rss_bytes'] is None or result['rss_samples'] < 1:
                raise RuntimeError(f'{name} has unavailable process-tree memory samples')
            report['stage_status'][name] = 'completed'
            return result
        stage('small_preflight', [sys.executable, ROOT/'scripts/reproduce_demo.py', '--binary', args.binary,
             '--expected-version', '0.5.0', '--label', 'source-preview', '--output', root/'preflight'])
        prep = module('case_prepare', ROOT/'benchmarks/evidence/prepare_representative.py')
        sources, _ = prep.source_lock(ROOT/protocol['source_lock'])
        cache = args.download_cache.resolve()
        report['cache_before'] = {row[0]: (cache/row[1]).exists() for row in sources}
        stage('acquisition_preparation', [sys.executable, ROOT/'benchmarks/evidence/prepare_representative.py',
              '--download-cache', cache, '--output', root/'prepared', '--contig', protocol['contig'],
              '--start', protocol['start'], '--end', protocol['end'], '--sites', protocol['probes']])
        prepared = read_json(root/'prepared/workloads.json')
        report['preparation'] = prepared['provenance']
        observed_preparation = {key: prepared['provenance'][key] for key in
                                ('window_records', 'window_max_read_length', 'dictionary_contigs')}
        if observed_preparation != {'window_records': 261761, 'window_max_read_length': 151,
                                    'dictionary_contigs': 195}:
            raise ValueError(f'pinned preparation observations changed: {observed_preparation}')
        report['preparation_observations_recomputed'] = observed_preparation
        report['acquired_bytes'] = sum(r['bytes'] for r in prepared['provenance']['sources']
                                        if not report['cache_before'][r['id']])
        report['authenticated_source_bytes'] = sum(r['bytes'] for r in prepared['provenance']['sources'])
        workload = next(row for row in prepared['workloads'] if row['id'] == 'window-bam-sites')
        for role in representatives().FILE_ROLES:
            workload[role]['path'] = str((root/'prepared'/workload[role]['path']).resolve())
        workload['cases'] = ['m256-t16384']
        workload['gates'] = {'min_admitted_budgets': 1, 'min_effective_tiles': 1, 'workers': [1]}
        workload['sample'] = 'HG002'
        import pysam
        with pysam.FastaFile(workload['reference']['path']) as fasta:
            contigs = [dict(name=n, length=l) for n,l in zip(fasta.references, fasta.lengths)]
        requests = freeze_queries(root/'prepared/window-sparse.vcf', root/'requests', contigs)
        report['stage_status']['query_freeze'] = 'completed'
        state = {'root': str(root), 'protocol': protocol, 'protocol_sha256': sha256(args.protocol),
                 'source_revision': args.expected_source, 'workload': workload, 'contigs': contigs, 'requests': requests,
                 'native': 'native/evidence-dataset.manifest.json', 'parquet': 'parquet', 'oracle': 'oracle.tsv',
                 'probe_sites': 'requests/probes.vcf', 'build_identity': 'build-identity.json',
                 'expected_native': {q['id']: 'native-results/'+q['id']+'.tsv' for q in requests}}
        write_json(root/'state.json', state)
        report['query_freeze'] = read_json(root/'requests/freeze.json')
        # Authenticate selected input identities before the oracle or timed extraction.
        inputs = {row['path']: {'sha256': row['sha256'], 'bytes': Path(row['path']).stat().st_size}
                  for row in (workload[role] for role in representatives().FILE_ROLES)}
        for path, record in inputs.items():
            if sha256(path) != record['sha256']:
                raise ValueError('prepared input differs from preparation identity')
        base = [sys.executable, Path(__file__).resolve(), '--binary', args.binary, '--state', root/'state.json']
        stage('oracle', [*base, '--worker', 'oracle', '--output', root/'oracle-work'])
        for request in requests:
            write_json(root/request['expected'], oracle_answers(root/'oracle.tsv', request))
        stage('materialize', [*base, '--worker', 'materialize', '--output', root/'materialize'])
        report['correctness'] = read_json(root/'materialize/result.json')
        stage('convert', [*base, '--worker', 'convert', '--output', root/'conversion'])
        state['parquet_inventory_sha256'] = sha256(root/'parquet/inventory.json')
        write_json(root/'state.json', state)
        report['stage_status']['trials'] = 'running'
        for repetition, order in enumerate(protocol['rounds'], 1):
            for path in order:
                destination = root/'trials'/f'r{repetition}-{path}'
                measured = measure([*base, '--worker', path, '--output', destination], root/'logs'/f'r{repetition}-{path}')
                measured.update(path=path, round=repetition, order=order,
                                observations=read_json(destination/'observations.json') if (destination/'observations.json').is_file() else [])
                report['trials'].append(measured)
                write_json(root/'report.json', report)
                if measured['exit_code'] or len(measured['prefixes']) != 3 or measured['rss_samples'] < 1:
                    raise RuntimeError(f'failed measured path r{repetition}-{path}')
        # Native A/B serialization must be equal at each fixed request, across every round.
        for request in requests:
            hashes = {sha256(root/'trials'/f'r{n}-{p}'/request['id']/'evidence.tsv') for n in (1,2,3) for p in ('A','B')}
            if len(hashes) != 1:
                raise ValueError('fresh/saved native TSV bytes differ')
        (root/'native-results').mkdir()
        for request in requests:
            source = root/'trials/r1-B'/request['id']/'evidence.tsv'
            destination = root/state['expected_native'][request['id']]
            shutil.copy2(source, destination)
            shutil.copy2(str(source)+'.manifest.json', str(destination)+'.manifest.json')
        report['synthetic_zero'] = synthetic_control(args.binary, root)
        report['stage_status']['trials'] = 'completed'
        report['stage_status']['controls'] = 'running'
        from controls import run_controls
        report['controls'] = run_controls(args.binary, root, root/'controls')
        report['stage_status']['controls'] = 'completed'
        report['stage_status']['source_rechecks'] = 'running'
        recheck_started = time.perf_counter()
        for path, record in inputs.items():
            if sha256(path) != record['sha256']:
                raise ValueError('immutable prepared input changed during study')
        for row in prepared['provenance']['sources']:
            if sha256(row['path']) != row['sha256']:
                raise ValueError('authenticated download changed during study')
        verify_inventory(root/'requests', report['query_freeze']['files'])
        authenticate(args.binary, args.build_identity, args.expected_source)
        report['source_rechecks'] = {'status': 'completed', 'seconds': time.perf_counter()-recheck_started,
                                    'prepared_inputs': inputs, 'binary_sha256': sha256(args.binary)}
        report['stage_status']['source_rechecks'] = 'completed'
        report['stage_status']['bundle'] = 'running'
        report['summary'] = summary(report)
        report['storage_logical_bytes'] = {'prepared_all_files': logical_bytes(root/'prepared'),
            'authenticated_downloads': report['authenticated_source_bytes'], 'native_minimum_replay': logical_bytes(root/'native'),
            'parquet_minimum_replay': logical_bytes(root/'parquet'), 'materialization_exports_and_receipts':
                sum(p.stat().st_size for p in (root/'materialize').iterdir() if p.is_file()),
            'all_query_outputs_and_logs': logical_bytes(root/'trials')}
        report['run_status'] = 'completed'
        report['complete_case_study_status'] = 'awaiting_clean_replay'
        write_json(root/'report.json', report)
        (root/'REPORT.md').write_text(human_report(report))
        with (root/'trials.csv').open('x') as stream:
            writer = csv.DictWriter(stream, fieldnames=['round','path','query','status','query_wall_seconds','verification_seconds','parent_prefix_seconds','answer_sha256'])
            writer.writeheader()
            for trial in report['trials']:
                for row, prefix in zip(trial['observations'], trial['prefixes']):
                    writer.writerow({key: (trial['round'] if key=='round' else prefix['parent_prefix_seconds'] if key=='parent_prefix_seconds' else row[key]) for key in writer.fieldnames})
        # Complete data bundle, without sources, duplicate caches or large build products.
        bundle = root/'replay-bundle'
        bundle.mkdir()
        for name in ('native', 'parquet', 'requests', 'expected', 'native-results', 'synthetic-zero'):
            shutil.copytree(root/name, bundle/name)
        for name in ('protocol-v1.json', 'build-identity.json', 'oracle.tsv', 'report.json', 'REPORT.md', 'trials.csv'):
            shutil.copy2(root/name, bundle/name)
        portable = {key: value for key,value in state.items() if key not in ('root','workload')}
        write_json(bundle/'state.json', portable)
        files = inventory(bundle)
        write_json(bundle/'inventory.json', {'files': files, 'source_revision': args.expected_source,
                   'protocol_sha256': sha256(args.protocol), 'build_identity_sha256': sha256(args.build_identity)})
        verify_inventory(bundle, files)
        retained = logical_bytes(bundle)
        if retained > protocol['bundle_limit_bytes']:
            raise ValueError('replay bundle exceeds fixed 128 MiB ceiling')
        report['stage_status']['bundle'] = 'completed'
        # Bundle's producer report describes scientific gates before packaging;
        # the outer artifact record authenticates the completed package.
        write_json(root/'report.json', report)
        write_json(root/'artifact.json', {'bundle_inventory_sha256': sha256(bundle/'inventory.json'), 'logical_bytes': retained,
                   'target_under_64_mib': retained < 64 << 20, 'source_revision': args.expected_source,
                   'binary_sha256': sha256(args.binary), 'retention_days': 30,
                   'trust': 'anchor this record to the exact hosted run; unsigned inventory alone does not authenticate authorship'})
    except Exception as error:
        report['run_status'] = 'failed'
        report['stage_status'] = {key: ('failed' if value == 'running' else value)
                                  for key,value in report['stage_status'].items()}
        report['failures'].append({'error': str(error), 'traceback': traceback.format_exc()})
        report['complete_case_study_status'] = 'failed_or_unstarted_gates'
        write_json(root/'report.json', report)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--build-identity', type=Path)
    parser.add_argument('--expected-source')
    parser.add_argument('--protocol', type=Path, default=HERE/'protocol-v1.json')
    parser.add_argument('--download-cache', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--worker', choices=('oracle','materialize','convert','A','B','D'))
    parser.add_argument('--state', type=Path)
    args = parser.parse_args()
    args.binary = args.binary.resolve(strict=True)
    args.output = args.output.resolve()
    if args.worker:
        if args.state is None:
            parser.error('--state required for worker')
        worker(args)
    else:
        if not all((args.build_identity, args.expected_source, args.download_cache)):
            parser.error('--build-identity, --expected-source and --download-cache are required')
        produce(args)


if __name__ == '__main__':
    main()
