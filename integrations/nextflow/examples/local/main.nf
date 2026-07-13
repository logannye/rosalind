nextflow.enable.dsl=2

include {
    ROSALIND_DOCTOR;
    ROSALIND_PLAN;
    ROSALIND_ANALYZE_SHARD;
    ROSALIND_MERGE;
    ROSALIND_VERIFY
} from '../../modules/rosalind/main'

workflow {
    shard_count = params.shard_count as int
    image = params.image
    meta = [id: params.id]
    indices = Channel.fromList((0..<shard_count).toList())

    doctor_input = Channel.value(tuple(meta, file(params.reference), file(params.bam), file(params.bai), image))
    ROSALIND_DOCTOR(doctor_input)

    plan_input = indices.map { index -> tuple(meta, file(params.reference), shard_count, index, params.budget_mb as int, image) }
    ROSALIND_PLAN(plan_input)

    analyses = ROSALIND_PLAN.out.plan.map { planned ->
        def (m, count, index, plan) = planned
        def prediction = new groovy.json.JsonSlurper().parseText(java.nio.file.Files.readString(plan))
        def plannedMiB = Math.ceil((prediction.predicted_peak_rss_bytes as long) / 1048576.0) as long
        tuple(m, file(params.reference), file(params.bam), file(params.bai), plan, plannedMiB,
              params.analyzer, params.format, count, index, params.budget_mb as int,
              params.scheduler_margin_mb as int, image)
    }
    ROSALIND_ANALYZE_SHARD(analyses)

    gathered = ROSALIND_ANALYZE_SHARD.out.artifacts.collect(flat: false).map { rows ->
        tuple(meta, rows.collect { it[2] }, rows.collect { it[3] }, image)
    }
    ROSALIND_MERGE(gathered)
    ROSALIND_VERIFY(ROSALIND_MERGE.out.merged.map {
        tuple(it[0], it[1], it[2], it[3] + it[4], image)
    })
}
