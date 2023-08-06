process FILTERINDELS {
    publishDir "${outdir}/3-variants/$pair_id", mode:'copy'

    input:
    tuple (val(pair_id), path(raw_indels))
    tuple val(name), path(ref), path(other_ref)
    val(round)
    val(outdir)

    output:
    tuple (val(pair_id),
	path("${pair_id}_filtered_indels_${round}.vcf"),
	path("${pair_id}_filtered_indels_${round}.vcf.idx"))

    script:
    """
    gatk VariantFiltration \
        -R $ref \
        -V $raw_indels \
        -O ${pair_id}_filtered_indels_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}
