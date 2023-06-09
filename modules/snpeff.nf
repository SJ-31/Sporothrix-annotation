process SNPEFF {
    publishDir "${outdir}/snpeff/$pair_id", mode:'copy'

    input:
    tuple (val(pair_id), path(filtered_snps),
    path(filtered_snps_index))
    val(outdir)

    output:
    path("*")

    script:
    """
    java -jar \$SNPEFF_JAR -v \
	-dataDir $params.snpeff_data \
	$params.snpeff_db \
	$filtered_snps > ${pair_id}_filtered_snps.ann.vcf
    """
}
