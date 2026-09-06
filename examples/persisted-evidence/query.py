#!/usr/bin/env python3
"""Verify a portable dataset, inspect a bounded query, and export Parquet parts."""
import argparse
from rosalind import open_dataset

parser = argparse.ArgumentParser()
parser.add_argument('manifest')
parser.add_argument('regions')
parser.add_argument('export_directory')
args = parser.parse_args()

dataset = open_dataset(args.manifest)
print(dataset.verify())
with dataset.batches(regions=args.regions, fields=['depths', 'alleles']) as run:
    rows = sum(batch.num_rows for batch in run)
    print({'positions': rows, 'receipt': str(run.result.manifest_path)})
exported = dataset.export_parquet(args.export_directory, regions=args.regions,
                                  fields=['depths', 'alleles', 'allele-quality'])
print({'parquet_directory': str(exported.directory), 'receipt': str(exported.manifest_path)})
