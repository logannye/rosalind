#!/usr/bin/env python3
"""Download a content-locked real-origin SAMtools example and prepare candidates.

Requires pysam==0.23.3 only for preparation. Candidate calls are illustrative,
not a benchmark truth set and not a Rosalind calling-accuracy measurement.
"""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            value.update(block)
    return value.hexdigest()


def main():
    import pysam
    import pysam.bcftools

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if pysam.__version__ != "0.23.3":
        parser.error("preparation is pinned to pysam==0.23.3")
    args.output.mkdir(parents=True, exist_ok=False)
    lock = json.loads(Path(__file__).with_name("sources.json").read_text())
    for resource in lock["files"]:
        destination = args.output / resource["name"]
        resource["url"] = ("https://raw.githubusercontent.com/samtools/samtools/"
                           f"{lock['source_commit']}/examples/{resource['name']}")
        with urllib.request.urlopen(resource["url"], timeout=60) as response:
            with destination.open("xb") as output:
                while block := response.read(65536):
                    output.write(block)
        if digest(destination) != resource["sha256"]:
            raise RuntimeError(f"source checksum mismatch: {destination}")
    reference = str(args.output / "ex1.fa")
    pysam.faidx(reference)
    raw_bam = str(args.output / "unsorted.bam")
    pysam.view("--no-PG", "-b", "-t", reference + ".fai", "-o", raw_bam,
               str(args.output / "ex1.sam.gz"), catch_stdout=False)
    bam = str(args.output / "sample.bam")
    pysam.sort("--no-PG", "-o", bam, raw_bam)
    pysam.index(bam)
    Path(raw_bam).unlink()
    likelihoods = str(args.output / "likelihoods.bcf")
    candidate_all = str(args.output / "candidate-all.vcf")
    mpileup = ["--no-version", "-f", reference, "-q", "20", "-Q", "20",
               "-Ob", "-o", likelihoods, bam]
    caller = ["--no-version", "-m", "-v", "--ploidy", "2", "-Ov", "-o",
              candidate_all, likelihoods]
    pysam.bcftools.mpileup(*mpileup, catch_stdout=False)
    pysam.bcftools.call(*caller, catch_stdout=False)
    candidates = args.output / "candidates.vcf"
    count = 0
    with open(candidate_all) as source, candidates.open("x") as output:
        for line in source:
            if line.startswith("#"):
                output.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields[3]) == 1 and fields[3] in "ACGT" and all(
                    len(alt) == 1 and alt in "ACGT" for alt in fields[4].split(",")):
                output.write(line)
                count += 1
    if count == 0:
        raise RuntimeError("independent caller produced no candidate SNVs")
    with pysam.FastaFile(reference) as fasta, (args.output / "targets.bed").open("x") as bed:
        for name, length in zip(fasta.references, fasta.lengths):
            bed.write(f"{name}\t0\t{length}\t{name}-full\n")
    lock.update({"pysam_version": pysam.__version__,
                 "htslib_samtools_bcftools_version": pysam.__samtools_version__,
                 "candidate_generator_argv": [["bcftools", "mpileup", *mpileup],
                                              ["bcftools", "call", *caller]],
                 "candidate_snvs": count,
                 "prepared_sha256": {path.name: digest(path) for path in args.output.iterdir()
                                     if path.is_file()}})
    (args.output / "preparation.json").write_text(json.dumps(lock, indent=2) + "\n")
    print(f"Prepared {count} independently generated candidate SNVs in {args.output}")


if __name__ == "__main__":
    main()
