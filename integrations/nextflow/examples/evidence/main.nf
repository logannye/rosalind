nextflow.enable.dsl=2

include { ROSALIND_EVIDENCE_PLAN; ROSALIND_PANEL_EVIDENCE; ROSALIND_VERIFY_EVIDENCE } from '../../modules/rosalind/evidence'

workflow {
    inputs = Channel.value(tuple([id: params.id], file(params.reference), file(params.bam),
                                 file(params.index), file(params.targets), params.budget_mb as int, params.image))
    ROSALIND_EVIDENCE_PLAN(inputs)
    analyses = ROSALIND_EVIDENCE_PLAN.out.planned.map { values ->
        def plan = new groovy.json.JsonSlurper().parseText(java.nio.file.Files.readString(values[7]))
        def mib = Math.ceil((plan.predicted_peak_rss_bytes as long) / 1048576.0) as long
        tuple(*values, mib, params.scheduler_margin_mb as int)
    }
    ROSALIND_PANEL_EVIDENCE(analyses)
    ROSALIND_VERIFY_EVIDENCE(ROSALIND_PANEL_EVIDENCE.out.artifacts)
}
