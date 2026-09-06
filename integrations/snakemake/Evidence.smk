import json
from pathlib import Path

SAMPLES = config["samples"]
BUDGET = int(config.get("budget_mb", 256))
MARGIN = int(config.get("scheduler_margin_mb", 64))
IMAGE = config["image"]

rule all:
    input:
        expand("results/{sample}/panel.tsv", sample=SAMPLES),
        expand("results/{sample}/positions.arrow", sample=SAMPLES),
        expand("results/{sample}/panel.tsv.manifest.json", sample=SAMPLES),
        expand("results/{sample}/evidence-verify.json", sample=SAMPLES)

rule evidence_plan:
    input:
        reference=lambda wc: SAMPLES[wc.sample]["reference"],
        bam=lambda wc: SAMPLES[wc.sample]["bam"],
        index=lambda wc: SAMPLES[wc.sample]["index"],
        targets=lambda wc: SAMPLES[wc.sample]["targets"]
    output: "work/{sample}/evidence-plan.json"
    log: "logs/{sample}/evidence-plan.log"
    params: budget=BUDGET
    container: IMAGE
    shell:
        "rosalind analyze panel-qc --reference-pack {input.reference:q} --alignments {input.bam:q} "
        "--alignment-index {input.index:q} --regions {input.targets:q} --memory-budget-mb {params.budget} "
        "--position-output planned-positions.arrow --output planned-panel.tsv --plan > {output:q} 2> {log:q}"

rule panel_evidence:
    input:
        reference=lambda wc: SAMPLES[wc.sample]["reference"],
        bam=lambda wc: SAMPLES[wc.sample]["bam"],
        index=lambda wc: SAMPLES[wc.sample]["index"],
        targets=lambda wc: SAMPLES[wc.sample]["targets"],
        plan="work/{sample}/evidence-plan.json"
    output:
        panel="results/{sample}/panel.tsv",
        positions="results/{sample}/positions.arrow",
        receipt="results/{sample}/panel.tsv.manifest.json"
    log: "logs/{sample}/panel-evidence.log"
    resources:
        mem_mb=lambda wildcards, input: (json.loads(Path(input.plan).read_text())["predicted_peak_rss_bytes"] + (1 << 20) - 1) // (1 << 20) + MARGIN
    params: budget=BUDGET
    container: IMAGE
    shell:
        "rosalind analyze panel-qc --reference-pack {input.reference:q} --alignments {input.bam:q} "
        "--alignment-index {input.index:q} --regions {input.targets:q} --min-callable-depth 10 "
        "--memory-budget-mb {params.budget} --position-output {output.positions:q} "
        "--output {output.panel:q} 2> {log:q}"

rule verify_evidence:
    input:
        panel="results/{sample}/panel.tsv",
        positions="results/{sample}/positions.arrow",
        receipt="results/{sample}/panel.tsv.manifest.json"
    output: "results/{sample}/evidence-verify.json"
    log: "logs/{sample}/evidence-verify.log"
    container: IMAGE
    shell: "rosalind verify --manifest {input.receipt:q} --json > {output:q} 2> {log:q}"
