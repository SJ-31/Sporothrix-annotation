process FILTERSNPS {
    publishDir "${outdir}/3-variants/$pair_id", mode:'copy'

    input:
    tuple val(pair_id), path(raw_snps)
    tuple val(name), path(ref), path(other_ref)
    val(round)
    val(outdir)

    output:
    tuple val(pair_id), path("${pair_id}_filtered_snps_${round}.vcf"), path("${pair_id}_filtered_snps_${round}.vcf.idx")

    script:
    """
    gatk VariantFiltration \
	-R $ref \
	-V $raw_snps \
	-O ${pair_id}_filtered_snps_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    """
}
