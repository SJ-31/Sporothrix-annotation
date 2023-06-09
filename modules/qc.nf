process QC {
    publishDir "$outdir/reports/$pair_id", mode: 'copy'

    input:
    tuple (val(pair_id),
        path("${pair_id}_dedup_metrics.txt"),
        path("${pair_id}_alignment_metrics.txt"),
        path("${pair_id}_insert_metrics.txt"),
        path("${pair_id}_insert_size_histogram.pdf"),
        path("${pair_id}_depth_out.txt"),
        path("${pair_id}_filtered_snps_1.vcf"),
        path("${pair_id}_filtered_snps_1.vcf.idx"),
        path("${pair_id}_filtered_snps_2.vcf"),
        path("${pair_id}_filtered_snps_2.vcf.idx"))
    val(outdir)

    output:
    path ("${pair_id}_report.csv")

    script:
    """
    parse_metrics.sh ${pair_id} > ${pair_id}_report.csv
    """
}
