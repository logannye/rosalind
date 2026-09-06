#!/usr/bin/env python3
"""Measure physical panel projection on a reproducible synthetic pressure case."""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys

from run import run_measured, sha256


def main():
    import pysam
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--budgets', default='32,64,96')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    budgets = [int(value) for value in args.budgets.split(',')]
    if args.repeats < 1 or any(value <= 0 for value in budgets):
        parser.error('positive budgets and repeats required')
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root / 'raw').mkdir()
    binary = args.binary.resolve()
    reference = root / 'reference.fa'
    length, step, read_length, copies = 200_000, 125, 250, 50
    reference.write_text('>chr1\n' + 'A' * length + '\n')
    pysam.faidx(str(reference))
    bam = root / 'panel.bam'
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': length}]}
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
        for start in range(0, length - read_length + 1, step):
            for copy in range(copies):
                read = pysam.AlignedSegment(out.header)
                read.query_name = f'r{start}-{copy}'
                read.query_sequence = ('C' if copy % 10 == 0 else 'A') * read_length
                read.query_qualities = [30 + copy % 10] * read_length
                read.flag = 16 if copy % 2 else 0
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigarstring = f'{read_length}M'
                out.write(read)
    pysam.index(str(bam))
    bed = root / 'panel.bed'
    bed.write_text(''.join(f'chr1\t{start}\t{start + 25000}\ttarget-{start}\n'
                           for start in range(0, length, 25000)))
    report = {'schema': 1, 'status': 'running', 'workload': 'synthetic-panel-projection-v1',
              'parameters': dict(length=length, step=step, read_length=read_length, copies=copies),
              'environment': dict(platform=platform.platform(), python=sys.version,
                  pysam=pysam.__version__, htslib=pysam.__samtools_version__,
                  binary_sha256=sha256(binary),
                  binary_version=subprocess.check_output([str(binary), '--version'], text=True).strip()),
              'inputs': {path.name: sha256(path) for path in
                  [reference, Path(str(reference) + '.fai'), bam, Path(str(bam) + '.bai'), bed]},
              'measurements': [],
              'limits': ['Synthetic pressure regression, not representative biological performance.',
                         'No OS-enforced limit is imposed.',
                         'Native extraction includes setup, hashing, analysis, encoding and publication; verification is separately timed.']}
    expected = None
    failed = False
    for repeat in range(args.repeats):
        # Alternate order to reduce systematic warm-cache advantage.
        for fields in (['all', 'depths,quality-sums'] if repeat % 2 == 0 else ['depths,quality-sums', 'all']):
            for budget in budgets:
                label = f'{fields.replace(",", "-")}-{budget}-{repeat + 1}'
                output = root / f'{label}.tsv'
                manifest = root / f'{label}.manifest.json'
                argv = [binary, 'analyze', 'panel-qc', '--reference', reference, '--alignments', bam,
                        '--regions', bed, '--fields', fields, '--memory-budget-mb', budget,
                        '--output', output, '--manifest', manifest]
                measured = run_measured(root, label, argv, root / f'{label}.stdout.txt')
                measured.update(fields=fields, budget_mib=budget, repeat=repeat + 1)
                if measured['exit_code'] == 0:
                    output_hash = sha256(output)
                    expected = output_hash if expected is None else expected
                    measured.update(output_bytes=output.stat().st_size, output_sha256=output_hash,
                                    output_matches=output_hash == expected,
                                    receipt_measurements=json.loads(manifest.read_text())['measurements'])
                    verified = run_measured(root, label + '-verify',
                        [binary, 'verify', '--manifest', manifest], root / f'{label}.verify.txt')
                    measured['verification'] = verified
                    failed |= output_hash != expected or verified['exit_code'] != 0
                else:
                    failed = True
                report['measurements'].append(measured)
                (root / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    report['status'] = 'failed' if failed else 'passed'
    (root / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    return int(failed)


if __name__ == '__main__':
    raise SystemExit(main())
