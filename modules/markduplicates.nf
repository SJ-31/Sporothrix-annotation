process MARKDUPLICATESSPARK {
    publishDir "${outdir}/dedup_sorted/$pair_id", mode:'copy'

    input:
    tuple val(pair_id), path(aligned_reads)
    val(outdir)

    output:
    tuple val(pair_id), path("${pair_id}_sorted_dedup.bam"), emit: sorted_bam
    tuple val(pair_id), path("${pair_id}_dedup_metrics.txt"), emit: dedup_qc_ch

    script:
    """
    mkdir -p ${params.tmpdir}/${workflow.runName}/${pair_id}
    gatk --java-options "-Djava.io.tmpdir=${params.tmpdir}/${workflow.runName}/${pair_id}" \
	MarkDuplicatesSpark \
	-I $aligned_reads \
	-M ${pair_id}_dedup_metrics.txt \
	-O ${pair_id}_sorted_dedup.bam
    rm -r ${params.tmpdir}/${workflow.runName}/${pair_id}
    """
}
