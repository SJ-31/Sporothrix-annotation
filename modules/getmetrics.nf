process GETMETRICS {
    publishDir "${outdir}/2-metrics/$pair_id", mode:'copy'

    input:
    tuple val(pair_id), path(sorted_dedup_reads)
    tuple val(name), path(ref), path(other_ref)
    val(outdir)

    output:
    tuple (val(pair_id),
	path("${pair_id}_alignment_metrics.txt"),
	path("${pair_id}_insert_metrics.txt"),
	path("${pair_id}_insert_size_histogram.pdf"),
	path("${pair_id}_depth_out.txt"))

    script:
    """
    java -jar \$PICARD_JAR \
        CollectAlignmentSummaryMetrics \
	R=${ref} \
        I=${sorted_dedup_reads} \
	O=${pair_id}_alignment_metrics.txt
    java -jar \$PICARD_JAR \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	OUTPUT=${pair_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${pair_id}_insert_size_histogram.pdf
    samtools depth -a ${sorted_dedup_reads} > ${pair_id}_depth_out.txt
    """
}
