nextflow.enable.dsl=2

process ROSALIND_EVIDENCE_PLAN {
    tag "${meta.id}"
    container { image }
    input:
    tuple val(meta), path(reference), path(bam), path(index), path(targets), val(budget_mb), val(image)
    output:
    tuple val(meta), path(reference), path(bam), path(index), path(targets), val(budget_mb), val(image), path('plan.json'), emit: planned
    script:
    """
    rosalind analyze panel-qc --reference-pack '${reference}' --alignments '${bam}' \
      --alignment-index '${index}' --regions '${targets}' \
      --memory-budget-mb ${budget_mb} --position-output planned-positions.arrow \
      --output planned-panel.tsv --plan > plan.json
    """
}

process ROSALIND_PANEL_EVIDENCE {
    tag "${meta.id}"
    container { image }
    memory { "${(planned_mb as long) + (margin_mb as long)} MB" }
    errorStrategy 'terminate'
    input:
    tuple val(meta), path(reference), path(bam), path(index), path(targets), val(budget_mb), val(image), path(plan), val(planned_mb), val(margin_mb)
    output:
    tuple val(meta), path('panel.tsv'), path('positions.arrow'), path('panel.tsv.manifest.json'), path(reference), path(bam), path(index), path(targets), val(image), emit: artifacts
    script:
    """
    rosalind analyze panel-qc --reference-pack '${reference}' --alignments '${bam}' \
      --alignment-index '${index}' --regions '${targets}' --min-callable-depth 10 \
      --memory-budget-mb ${budget_mb} --position-output positions.arrow --output panel.tsv
    """
}

process ROSALIND_VERIFY_EVIDENCE {
    tag "${meta.id}"
    container { image }
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(meta), path(panel), path(positions), path(receipt), path(reference), path(bam), path(index), path(targets), val(image)
    output:
    tuple val(meta), path('panel.tsv'), path('positions.arrow'), path('panel.tsv.manifest.json'), path('verify.json'), emit: verified
    script:
    """
    rosalind verify --manifest '${receipt}' --json > verify.json
    """
}
