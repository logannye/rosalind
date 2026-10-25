"""Deterministic saved-evidence refusal controls, never timed as queries."""
from __future__ import annotations

from pathlib import Path
import shutil

from case_common import (command, inventory, native_query, native_table, probes,
                         read_json, safe_path, sha256, write_json, write_sites)


def publication_disposition(directory, output):
    """A nonzero process cannot publish a successful primary output or receipt."""
    directory, output = Path(directory), Path(output)
    files = inventory(directory)
    receipt = Path(str(output) + '.manifest.json')
    receipt_status = None
    if receipt.exists():
        value = read_json(receipt)
        receipt_status = value.get('params', {}).get('run_status')
        if receipt_status not in ('failed', 'refused', 'resource-failed', 'interrupted'):
            raise ValueError('negative control published an unclassified or successful receipt')
    if output.exists():
        raise ValueError('negative control published a final primary output')
    return {'final_primary_exists': False, 'receipt_status': receipt_status,
            'files': files, 'success_published': False}


def _negative(binary, manifest, selection, directory, *, regions=False):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    output = directory / 'attempt.tsv'
    argv = native_query(binary, manifest, selection, output)
    if regions:
        argv[argv.index('--sites')] = '--regions'
    result = command(argv, directory / 'command', expect_failure=True)
    result.update(publication_disposition(directory, output))
    result['diagnostic'] = (directory / 'command/stderr.log').read_text(errors='replace')[-65536:]
    return result


def consumed_partition(bundle, state):
    """Select the first Q1 locus's ownership partition, using no read counts."""
    manifest = safe_path(bundle, state['native'])
    descriptor = read_json(manifest.parent / 'dataset.descriptor.json')
    query = next(request for request in state['requests'] if request['id'] == 'Q1')
    locus = query['rows'][0]
    ranks = {item['name']: index for index, item in enumerate(state['contigs'])}
    rank, position = ranks[locus['contig']], locus['pos'] - 1
    matches = [part for part in descriptor['partitions'] if part['contig'] == rank
               and any(interval['start'] <= position < interval['end']
                       for interval in part['intervals'])]
    if len(matches) != 1:
        raise ValueError('Q1 locus does not identify exactly one consumed partition')
    payload = safe_path(manifest.parent, matches[0]['arrow']['path'])
    if not payload.is_file() or payload.stat().st_size < 1:
        raise ValueError('consumed native evidence payload is absent or empty')
    return payload, locus


def run_controls(binary, bundle, output):
    """Keep the intact bundle immutable; mutate only a newly owned complete copy."""
    bundle, output = Path(bundle), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    report = {'schema_version': 1, 'status': 'running', 'checks': {}}
    try:
        state = read_json(bundle / 'state.json')
        manifest = safe_path(bundle, state['native'])
        baseline = inventory(manifest.parent)
        selections_before = inventory(safe_path(bundle, 'requests'))
        selected = probes(safe_path(bundle, state['probe_sites']))
        if not selected:
            raise ValueError('missing frozen probes')
        mask = {(row['contig'], row['pos']) for row in selected}
        protocol = read_json(bundle / 'protocol-v1.json')
        absent = next(pos for pos in range(protocol['start'] + 1, protocol['end'] + 1)
                      if (protocol['contig'], pos) not in mask)
        absent_file = output / 'absent.bed'
        absent_file.write_text(f"{protocol['contig']}\t{absent-1}\t{absent}\n")
        expected_absent = safe_path(bundle, 'requests/absent.bed')
        if absent_file.read_bytes() != expected_absent.read_bytes():
            raise ValueError('frozen absent-locus control differs from probe mask')
        report['checks']['missing_locus'] = _negative(binary, manifest, absent_file,
                                                     output / 'missing', regions=True)
        if 'absent from the dataset' not in report['checks']['missing_locus']['diagnostic']:
            raise ValueError('missing-locus control did not reach the coverage refusal')
        wrong = dict(selected[0])
        wrong['ref'] = next(base for base in 'ACGT' if base not in (wrong['ref'], wrong['alt']))
        wrong_file = output / 'wrong-ref.vcf'
        write_sites(wrong_file, [wrong], state['contigs'])
        if wrong_file.read_bytes() != safe_path(bundle, 'requests/wrong-ref.vcf').read_bytes():
            raise ValueError('frozen REF control differs from probe mask')
        report['checks']['wrong_ref'] = _negative(binary, manifest, wrong_file, output / 'wrong-ref')
        if 'query SNV REF differs' not in report['checks']['wrong_ref']['diagnostic']:
            raise ValueError('REF control did not reach consumed-row REF validation')

        import checked_parquet
        parquet = safe_path(bundle, state['parquet'])
        for name, row, diagnostic in (
                ('parquet_missing', dict(selected[0], contig=protocol['contig'], pos=absent), 'missing requested locus'),
                ('parquet_wrong_ref', wrong, 'REF mismatch')):
            try:
                checked_parquet.query(parquet, [row], expected_inventory_sha256=state['parquet_inventory_sha256'])
            except ValueError as error:
                if diagnostic not in str(error):
                    raise ValueError(f'{name} failed for an unrelated reason') from error
                report['checks'][name] = {'status': 'expected_negative_passed', 'diagnostic': str(error),
                                        'success_published': False}
            else:
                raise ValueError(f'{name} returned a successful result')

        source, locus = consumed_partition(bundle, state)
        corrupted = output / 'corrupted-native'
        shutil.copytree(manifest.parent, corrupted, symlinks=False)
        target = corrupted / source.relative_to(manifest.parent)
        original_hash = sha256(target)
        offset = target.stat().st_size // 2
        with target.open('r+b') as payload:
            payload.seek(offset)
            original_byte = payload.read(1)[0]
            payload.seek(offset)
            payload.write(bytes([original_byte ^ 1]))
        changed_hash = sha256(target)
        if changed_hash == original_hash:
            raise ValueError('fault injection did not change the consumed payload')
        changed = [name for name, item in inventory(corrupted).items() if baseline.get(name) != item]
        if changed != [source.relative_to(manifest.parent).as_posix()]:
            raise ValueError('fault control changed more than the selected evidence payload')
        bad_manifest = corrupted / manifest.name
        verified = command([binary, 'dataset', 'verify', '--dataset', bad_manifest,
                            '--memory-budget-mb', '256'], output / 'tamper-verify', expect_failure=True)
        q1 = next(request for request in state['requests'] if request['id'] == 'Q1')
        refused = _negative(binary, bad_manifest, safe_path(bundle, q1['sites']), output / 'tamper-query')
        report['checks']['tamper'] = {'status': 'expected_negative_passed',
            'consumed_locus': locus, 'payload': source.relative_to(manifest.parent).as_posix(),
            'offset': offset, 'original_byte': original_byte, 'changed_byte': original_byte ^ 1,
            'before_sha256': original_hash, 'after_sha256': changed_hash,
            'manifests_unchanged': True, 'dataset_verify': verified, 'query': refused,
            'verification_diagnostic': (output / 'tamper-verify/stderr.log').read_text(errors='replace')[-65536:]}
        report['checks']['synthetic_zero'] = synthetic_zero(binary, bundle, output / 'synthetic-zero')
        if inventory(manifest.parent) != baseline or inventory(safe_path(bundle, 'requests')) != selections_before:
            raise ValueError('intact baseline evidence or requests changed during negative controls')
        report.update(status='expected_negative_passed', intact_native_and_requests_unchanged=True)
    except BaseException as error:
        report.update(status='failed', error=f'{type(error).__name__}: {error}')
        write_json(output / 'report.json', report)
        raise
    write_json(output / 'report.json', report)
    return report


def synthetic_zero(binary, bundle, output):
    """Replay the existing authored specimen-C zero separately from HG002."""
    root = safe_path(bundle, 'synthetic-zero')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    if ((root / 'zero.bed').read_text() != 'synthetic\t39\t40\n'
            or (root / 'absent.bed').read_text() != 'synthetic\t40\t41\n'):
        raise ValueError('synthetic zero/missing controls differ from authored coordinates')
    expected_path = Path(__file__).resolve().parents[2] / 'examples/cohort-reanalysis/expected.json'
    if sha256(root / 'authored-expected.json') != sha256(expected_path):
        raise ValueError('authored synthetic expectations differ from exact source')
    authored = read_json(root / 'authored-expected.json')
    expected_row = next(row for row in authored['members']['specimen-c']['full_rows'] if row[0] == 40)
    manifest = root / 'native/evidence-dataset.manifest.json'
    target = output / 'zero.tsv'
    argv = native_query(binary, manifest, root / 'zero.bed', target)
    argv[argv.index('--sites')] = '--regions'
    extraction = command(argv, output / 'zero-query')
    verification = command([binary, 'verify', '--manifest', str(target) + '.manifest.json', '--json'],
                           output / 'zero-verify')
    rows = list(native_table(target))
    if (len(rows) != 1 or rows[0]['contig'] != 'synthetic' or rows[0]['ref'] != 'A'
            or [int(rows[0][column]) for column in authored['full_row_columns']] != expected_row
            or any(int(rows[0][base + strand]) != 0 for base in 'acgt' for strand in ('_fwd', '_rev'))):
        raise ValueError('stored synthetic zero differs from hand-authored expectation')
    missing = _negative(binary, manifest, root / 'absent.bed', output / 'unmeasured', regions=True)
    if 'absent from the dataset' not in missing['diagnostic']:
        raise ValueError('synthetic missing control did not reach coverage refusal')
    return {'status': 'completed', 'label': 'authored synthetic fixture, separate from HG002',
            'stored_zero_position': 40, 'unmeasured_position': 41, 'callable_depth': 0,
            'query': extraction, 'verification': verification, 'unmeasured': missing}
