nextflow.enable.dsl=2

process ROSALIND_DOCTOR {
    tag "${meta.id}"
    container { image }
    input:
    tuple val(meta), path(reference), path(bam), path(bai), val(image)
    output:
    tuple val(meta), path('doctor.json'), emit: report
    script:
    """
    rosalind doctor --reference-pack ${reference} --alignments ${bam} --json > doctor.json
    """
}

process ROSALIND_PLAN {
    tag "${meta.id}:shard-${shard_index}"
    container { image }
    input:
    tuple val(meta), path(reference), val(shard_count), val(shard_index), val(budget_mb), val(image)
    output:
    tuple val(meta), val(shard_count), val(shard_index), path('plan.json'), emit: plan
    script:
    """
    rosalind plan --reference-pack ${reference} --shard-count ${shard_count} \
      --shard-index ${shard_index} --budget-mb ${budget_mb} --json > plan.json
    """
}

process ROSALIND_ANALYZE_SHARD {
    tag "${meta.id}:shard-${shard_index}"
    container { image }
    errorStrategy 'terminate'
    memory { "${(planned_memory_mb as long) + (scheduler_margin_mb as long)} MB" }
    input:
    tuple val(meta), path(reference), path(bam), path(bai), path(plan), val(planned_memory_mb),
          val(analyzer), val(format), val(shard_count), val(shard_index), val(budget_mb),
          val(scheduler_margin_mb), val(image)
    output:
    tuple val(meta), val(shard_index), path("shard-${shard_index}.${format == 'arrow-ipc' ? 'arrow' : 'tsv'}"),
          path("shard-${shard_index}.${format == 'arrow-ipc' ? 'arrow' : 'tsv'}.manifest.json"), emit: artifacts
    script:
    def output = "shard-${shard_index}.${format == 'arrow-ipc' ? 'arrow' : 'tsv'}"
    def command = analyzer == 'features' ? "features --format ${format}" : "analyze ${analyzer}"
    """
    rosalind ${command} --reference-pack ${reference} --alignments ${bam} \
      --shard-count ${shard_count} --shard-index ${shard_index} \
      --memory-budget-mb ${budget_mb} --enforce --output ${output}
    """
}

process ROSALIND_MERGE {
    tag "${meta.id}"
    container { image }
    input:
    tuple val(meta), path(artifacts), path(manifests), val(image)
    output:
    tuple val(meta), path('merged.arrow'), path('merged.arrow.manifest.json'),
          path(artifacts), path(manifests), emit: merged
    script:
    """
    args=()
    for receipt in ${manifests}; do args+=(--manifest "\$receipt"); done
    rosalind merge "\${args[@]}" --inputs . --output merged.arrow
    """
}

process ROSALIND_VERIFY {
    tag "${meta.id}"
    container { image }
    input:
    tuple val(meta), path(artifact), path(manifest), path(evidence_inputs), val(image)
    output:
    tuple val(meta), path(artifact), path(manifest), path('verify.json'), emit: verified
    script:
    """
    rosalind verify --manifest ${manifest} --json > verify.json
    """
}

process ROSALIND_REPRODUCE {
    tag "${meta.id}"
    container { image }
    input:
    tuple val(meta), path(artifact), path(manifest), path(inputs), val(image)
    output:
    tuple val(meta), path(artifact), path(manifest), path('reproduce.json'), path('*.repro.json'), emit: reproduced
    script:
    """
    rosalind reproduce --manifest ${manifest} --inputs ${inputs} --json > reproduce.json
    """
}
