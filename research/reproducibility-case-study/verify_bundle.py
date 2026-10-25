#!/usr/bin/env python3
"""Recompute a portable case study using saved evidence and an exact binary only."""
from __future__ import annotations

import argparse
from pathlib import Path
import re

import build_identity
import checked_parquet
from case_common import (command, freeze_queries, inventory, load_protocol, native_query, native_table,
                         oracle_answers, read_json, reduce_native, safe_path, sha256,
                         verify_inventory, write_json)
from controls import run_controls


def parquet_oracle_comparison(directory, oracle, expected_inventory_sha256):
    """Recheck all 34 physical fields against raw-oracle observations in batches."""
    import pyarrow.parquet as pq
    directory = Path(directory)
    before = checked_parquet.verify(directory, expected_inventory_sha256=expected_inventory_sha256)
    table = read_json(directory / 'inventory.json')
    rows = iter(native_table(oracle))
    count = 0
    for batch in pq.ParquetFile(safe_path(directory, table['files'][0]['path'])).iter_batches(
            batch_size=1024, use_threads=False):
        for physical in batch.to_pylist():
            try:
                expected = next(rows)
            except StopIteration as error:
                raise ValueError('Parquet has extra rows beyond independent oracle') from error
            converted = {'contig': expected['contig'], 'pos': int(expected['pos']),
                         'ref': expected['ref'], 'requested_alts': expected['requested_alts'].replace(',', '').replace('.', '')}
            converted.update({key: int(expected[key]) for key in checked_parquet.SCALARS})
            for name, size in checked_parquet.HISTOGRAMS.items():
                bins = [0] * size
                if expected[name] != '.':
                    previous = -1
                    for item in expected[name].split(','):
                        index, value = map(int, item.split(':'))
                        if not previous < index < size or not 1 <= value <= 2**64-1:
                            raise ValueError('invalid independent oracle histogram')
                        bins[index] = value
                        previous = index
                converted[name] = bins
            if physical != converted:
                raise ValueError(f'Parquet legacy fields differ from oracle at row {count+1}')
            count += 1
    if next(rows, None) is not None:
        raise ValueError('Parquet omits independent oracle rows')
    checked_parquet.verify(directory, expected_inventory_sha256=expected_inventory_sha256)
    return {'rows': count, 'columns': 34, 'mismatches': 0, 'full_schema_values_equal': True,
            'inventory_sha256': before['inventory_sha256']}


def inspect_bundle(bundle, expected_inventory_sha256=None):
    """Validate completeness before trusting any state or producer pass field."""
    bundle = Path(bundle)
    if bundle.is_symlink() or not bundle.is_dir():
        raise ValueError('bundle must be a regular directory')
    bundle = bundle.resolve(strict=True)
    path = bundle / 'inventory.json'
    if path.is_symlink():
        raise ValueError('inventory symlink refused')
    digest = sha256(path)
    if expected_inventory_sha256 is not None:
        if not re.fullmatch('[0-9a-f]{64}', expected_inventory_sha256) or digest != expected_inventory_sha256:
            raise ValueError('bundle inventory differs from retained external digest')
    record = read_json(path)
    files = record['files']
    verify_inventory(bundle, files)
    if inventory(bundle, exclude=('inventory.json',)) != files:
        raise ValueError('bundle file set differs from its complete inventory')
    if sum(item['bytes'] for item in files.values()) > 128 << 20:
        raise ValueError('bundle exceeds 128 MiB envelope')
    state = read_json(bundle / 'state.json')
    protocol = load_protocol(bundle / 'protocol-v1.json')
    if (state['protocol_sha256'] != sha256(bundle / 'protocol-v1.json')
            or record['protocol_sha256'] != state['protocol_sha256']
            or record['source_revision'] != state['source_revision']
            or record['build_identity_sha256'] != sha256(safe_path(bundle, state['build_identity']))):
        raise ValueError('bundle source/protocol/build binding disagrees')
    required = [state['native'], state['oracle'], state['probe_sites'], state['build_identity'],
                'requests/freeze.json', 'requests/absent.bed', 'requests/wrong-ref.vcf',
                state['parquet'] + '/inventory.json', 'synthetic-zero/native/evidence-dataset.manifest.json',
                'synthetic-zero/zero.bed', 'synthetic-zero/absent.bed', 'synthetic-zero/authored-expected.json',
                'report.json', 'REPORT.md', 'trials.csv']
    if [request['id'] for request in state['requests']] != ['Q1', 'Q2', 'Q3']:
        raise ValueError('bundle requires exactly Q1/Q2/Q3 in order')
    for request in state['requests']:
        expected_native = state['expected_native'][request['id']]
        required += [request['sites'], request['expected'], expected_native,
                     expected_native + '.manifest.json']
    for name in required:
        if name not in files or not safe_path(bundle, name).is_file():
            raise ValueError(f'missing required bundle artifact: {name}')
    for leaf in (safe_path(bundle, state['native']).parent, bundle / 'synthetic-zero/native'):
        descriptor_path = leaf / 'dataset.descriptor.json'
        descriptor = read_json(descriptor_path)
        for partition in descriptor['partitions']:
            for role in ('arrow', 'receipt'):
                path = safe_path(leaf, partition[role]['path'])
                if path.relative_to(bundle).as_posix() not in files or not path.is_file():
                    raise ValueError('bundle omits a descriptor-required partition artifact')
    if sha256(safe_path(bundle, state['parquet'] + '/inventory.json')) != state['parquet_inventory_sha256']:
        raise ValueError('Parquet inventory differs from outer run binding')
    for name in files:
        if name.lower().endswith(('.bam', '.bai', '.cram', '.crai', '.csi', '.fa', '.fasta', '.fai', '.fa.gz')):
            raise ValueError('saved-only bundle includes a raw alignment/reference input')
    return state, record, protocol, digest


def reconstruct_requests(bundle, state, output):
    """Reconstruct requests from retained probes; never trust changed definitions."""
    generated = freeze_queries(safe_path(bundle, state['probe_sites']), output, state['contigs'])
    if generated != state['requests']:
        raise ValueError('retained requests differ from independently reconstructed frozen definitions')
    for path in Path(output).iterdir():
        if path.name == 'freeze.json':
            continue
        retained = safe_path(bundle, 'requests/' + path.name)
        if path.read_bytes() != retained.read_bytes():
            raise ValueError(f'retained request/control bytes changed: {path.name}')
    freeze = read_json(safe_path(bundle, 'requests/freeze.json'))
    if freeze.get('requests') != generated:
        raise ValueError('request freeze record disagrees with regenerated requests')
    verify_inventory(safe_path(bundle, 'requests'), freeze['files'])
    answers = {}
    for request in generated:
        expected = oracle_answers(safe_path(bundle, state['oracle']), request)
        if read_json(safe_path(bundle, request['expected'])) != expected:
            raise ValueError('retained expected answers differ from independent oracle derivation')
        answers[request['id']] = expected
    return generated, answers


def source_absence(bundle, state):
    manifest = safe_path(bundle, state['native'])
    descriptor = read_json(manifest.parent / 'dataset.descriptor.json')
    checks = []
    for source in descriptor['sources']:
        if source['role'] in ('alignments', 'alignment-index', 'reference', 'reference-fai',
                              'cram-reference', 'cram-reference-fai'):
            path = Path(source['path'])
            checks.append({'role': source['role'], 'original_path': str(path), 'exists': path.exists()})
    original = read_json(safe_path(bundle, state['expected_native']['Q1'] + '.manifest.json'))
    for source in original.get('inputs', []):
        path = Path(source['path'])
        checks.append({'role': 'original-saved-query-input', 'original_path': str(path), 'exists': path.exists()})
    if not checks or any(check['exists'] for check in checks):
        raise ValueError('saved-only replay can still access an original raw/cache/query input')
    return {'checks': checks, 'original_sources_and_old_cache_absent': True,
            'scope': 'observed absent original paths plus supplied bundle inventory; separate clean CI job is the isolation boundary'}


def verify_bundle(bundle, binary, identity_path, output, expected_inventory_sha256=None):
    bundle, binary, identity_path, output = map(Path, (bundle, binary, identity_path, output))
    if output.resolve().is_relative_to(bundle.resolve()):
        raise ValueError('replay output must be outside the immutable bundle')
    output.mkdir(parents=True, exist_ok=False)
    report = {'schema_version': 1, 'status': 'running', 'commands': [], 'queries': [],
              'trust_boundary': 'consistency with retained source/build/inventory; coordinated replacement of all provenance is outside scope'}
    write_json(output / 'report.json', report)
    try:
        state, bundle_inventory, protocol, initial_digest = inspect_bundle(bundle, expected_inventory_sha256)
        report.update(source_revision=state['source_revision'], protocol_sha256=state['protocol_sha256'],
                      bundle_inventory_sha256=initial_digest,
                      externally_anchored_inventory=expected_inventory_sha256 is not None)
        if sha256(identity_path) != sha256(safe_path(bundle, state['build_identity'])):
            raise ValueError('selected external build identity differs from bundled identity')
        identity = build_identity.authenticate(binary, identity_path, state['source_revision'])
        report['build_identity'] = identity
        report['verifier_environment'] = build_identity.environment_identity()
        report['source_absence'] = source_absence(bundle, state)
        requests, expected = reconstruct_requests(bundle, state, output / 'reconstructed-requests')
        manifest = safe_path(bundle, state['native'])
        result = command([binary, 'dataset', 'verify', '--dataset', manifest, '--memory-budget-mb', '256'],
                         output / 'native-full-verification')
        report['commands'].append(result)
        native_verified = read_json(output / 'native-full-verification/stdout.log')
        if native_verified.get('verified_loci') != protocol['probes'] or native_verified.get('original_sources_rehashed') is not False:
            raise ValueError('full native dataset verification counters disagree')
        report['native_verification'] = native_verified
        full_output = output / 'full-legacy.tsv'
        full_argv = native_query(binary, manifest, safe_path(bundle, state['probe_sites']), full_output)
        full_argv[full_argv.index('--fields') + 1] = 'all'
        report['commands'].append(command(full_argv, output / 'full-legacy-extract'))
        report['commands'].append(command([binary, 'verify', '--manifest', str(full_output) + '.manifest.json',
                                           '--json'], output / 'full-legacy-verify'))
        if full_output.read_bytes() != safe_path(bundle, state['oracle']).read_bytes():
            raise ValueError('complete legacy saved evidence differs from independent raw oracle')
        report['full_legacy_comparison'] = {'rows': protocol['probes'], 'columns': 34, 'mismatches': 0,
                                            'exact_tsv_bytes_equal': True}
        report['parquet_verification'] = checked_parquet.verify(safe_path(bundle, state['parquet']),
            expected_inventory_sha256=state['parquet_inventory_sha256'])
        report['parquet_oracle_comparison'] = parquet_oracle_comparison(safe_path(bundle, state['parquet']),
            safe_path(bundle, state['oracle']), state['parquet_inventory_sha256'])
        for request in requests:
            name = request['id']
            native_output = output / (name + '.tsv')
            launched = command(native_query(binary, manifest, safe_path(bundle, request['sites']), native_output),
                               output / (name + '-native'))
            verified = command([binary, 'verify', '--manifest', str(native_output) + '.manifest.json', '--json'],
                               output / (name + '-verify'))
            report['commands'].extend([launched, verified])
            receipt = read_json(Path(str(native_output) + '.manifest.json'))
            values = receipt['measurements']
            if (receipt.get('params', {}).get('run_status') != 'completed'
                    or values.get('execution.alignment_record_visits') != '0'
                    or values.get('original_sources_rehashed') != 'false'):
                raise ValueError('saved query lacks completed zero-alignment-visit provenance')
            native = reduce_native(native_output, request)
            if native != expected[name]:
                raise ValueError(f'{name}: native saved answer differs from independent oracle')
            if native_output.read_bytes() != safe_path(bundle, state['expected_native'][name]).read_bytes():
                raise ValueError(f'{name}: native output bytes changed after relocation')
            parquet, metrics = checked_parquet.query(safe_path(bundle, state['parquet']), request['rows'],
                thresholds=request['thresholds'], expected_inventory_sha256=state['parquet_inventory_sha256'])
            if parquet != expected[name]:
                raise ValueError(f'{name}: checked Parquet answer differs from independent oracle')
            write_json(output / (name + '-canonical.json'), native)
            report['queries'].append({'id': name, 'rows': len(native), 'canonical_equal': True,
                'native_bytes_equal': True, 'native_sha256': sha256(native_output), 'parquet': metrics,
                'alignment_record_visits': 0, 'original_sources_rehashed': False})
            write_json(output / 'report.json', report)
        q1_receipt = safe_path(bundle, state['expected_native']['Q1'] + '.manifest.json')
        report['commands'].append(command([binary, 'reproduce', '--manifest', q1_receipt,
            '--inputs', bundle, '--binary', binary, '--json', '--no-attest'], output / 'native-reproduce'))
        replay = read_json(output / 'native-reproduce/stdout.log')
        if replay.get('verdict') != 'REPRODUCED' or replay.get('code_identity_matches') is not True:
            raise ValueError('original native Q1 receipt did not reproduce with matching code')
        report['native_reproduction'] = replay
        report['controls'] = run_controls(binary, bundle, output / 'negative-controls')
        final_state, _, _, final_digest = inspect_bundle(bundle, initial_digest)
        if final_state != state or sha256(binary) != identity['binary_sha256']:
            raise ValueError('bundle or executable changed during replay')
        report.update(status='completed', scientific_comparisons_recomputed=True,
                      bundle_unchanged=final_digest == initial_digest,
                      supplied_files=bundle_inventory['files'])
    except BaseException as error:
        report.update(status='failed', error=f'{type(error).__name__}: {error}')
        write_json(output / 'report.json', report)
        raise
    write_json(output / 'report.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--build-identity', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-inventory-sha256')
    args = parser.parse_args()
    verify_bundle(args.bundle, args.binary, args.build_identity, args.output, args.expected_inventory_sha256)


if __name__ == '__main__':
    main()
