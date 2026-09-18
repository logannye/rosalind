#!/usr/bin/env python3
"""Render the retained baseline measurements; requires matplotlib==3.10.6."""
from pathlib import Path
import csv
import json
import statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
report = json.loads((ROOT / 'report.json').read_text())
audit = json.loads((ROOT / 'audit.json').read_text())
measured = {m['label']: m for m in audit['measurements']}
verification = {(x['workload'], x['case'], x['phase']): x for x in audit['verification_summaries']}
with (ROOT / 'medians.csv').open('w') as output:
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(['workload', 'case', 'phase', 'repeats', 'wall_seconds', 'rss_mib', 'verify_seconds'])
    for row in sorted(report['summaries'], key=lambda x: (x['workload'], x['case'], x['phase'])):
        key = (row['workload'], row['case'], row['phase'])
        check = verification.get(key)
        writer.writerow([*key, row['completed_repeats'], row['wall_seconds']['median'],
                         row['peak_rss_bytes']['median'] / 2**20,
                         check['wall_seconds']['median'] if check else ''])

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(2, 2, figsize=(10.5, 7), constrained_layout=True)
cases = ['m96-t128', 'm128-t4096', 'm256-t16384']
for row_index, (selection, title) in enumerate([('sites', '1,000 sparse loci'), ('regions', '10,000 panel positions')]):
    for encoding, color in [('bam', '#2563a5'), ('cram', '#bc5c16')]:
        workload = f'window-{encoding}-{selection}'
        data = {s['case']: s for s in report['summaries'] if s['workload'] == workload and s['phase'] == 'fresh'}
        for column, (metric, divisor, ylabel) in enumerate([
            ('wall_seconds', 1, 'End-to-end wall time (s)'), ('peak_rss_bytes', 2**20, 'Process peak RSS (MiB)')]):
            ax = axes[row_index, column]
            values = [data[c][metric]['median'] / divisor for c in cases]
            errors = [[(data[c][metric]['median'] - data[c][metric]['minimum']) / divisor for c in cases],
                      [(data[c][metric]['maximum'] - data[c][metric]['median']) / divisor for c in cases]]
            ax.errorbar(range(3), values, yerr=errors, label=encoding.upper(), color=color,
                        marker='o', linewidth=1.8, capsize=3)
            if column == 1:
                predictions = [statistics.median(measured[m['label']]['predicted_peak_rss_bytes']
                                for m in report['measurements'] if m['workload'] == workload
                                and m.get('case') == case and m['phase'] == 'fresh') / divisor for case in cases]
                ax.plot(range(3), predictions, color=color, linestyle=':', alpha=.75,
                        label=encoding.upper() + ' plan')
            ax.set_title(title)
            ax.set_ylabel(ylabel)
            ax.set_xticks(range(3), ['96 / 128', '128 / 4,096', '256 / 16,384'])
            ax.set_xlabel('Declared budget (MiB) / requested tile (bases)')
            ax.grid(axis='y', alpha=.18)
    for ax in axes[row_index]:
        ax.set_ylim(bottom=0)
    axes[row_index, 0].legend(frameon=False)
    axes[row_index, 1].legend(frameon=False, ncols=2)
fig.suptitle('HG002 baseline: budget and tile configuration sweep\nThree runs per point; bars show min–max; one worker', fontsize=14)
fig.savefig(ROOT / 'resource-curve.png', dpi=160)
fig.savefig(ROOT / 'resource-curve.svg')

# Matplotlib SVG path data may end in spaces; normalize formatting only.
svg = ROOT / "resource-curve.svg"
svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
