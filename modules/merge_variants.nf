process MERGE_VARIANTS {
    publishDir "$outdir/variants/$pair_id", mode: 'copy'
    input:
    tuple (val(pair_id),
	path(snps),
	path(snps_index))
    tuple (val(pair_id),
	path(indels),
	path(indel_index))
    val(outdir)
    //
    output:
    tuple (val(pair_id), path("${pair_id}_merged.vcf.gz"))
    //
    script:
    """
    bgzip $snps
    bgzip $indels
    bcftools index ${snps}.gz
    bcftools index ${indels}.gz
    bcftools concat -a ${snps}.gz ${indels}.gz -O b \
    > ${pair_id}_merged.vcf.gz
    """
    //
}
