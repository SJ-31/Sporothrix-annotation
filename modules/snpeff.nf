process SNPEFF {
    publishDir "${outdir}/$pair_id", mode:'copy'

    input:
    tuple (val(pair_id), path(merged_variants),
    path(filtered_snps_index))
    val(outdir)

    output:
    path("${pair_id}_snpeEff*")
    path("${pair_id}_annotated_vars.vcf")

    script:
    """
    snpEff_wrapper.sh $merged_variants \
    $params.snpeff_database \
    $params.snpeff_db_name \
    $pair_id
    mv snpEff_genes.txt ${pair_id}_snpEff_genes.txt
    mv snpEff_summary.html ${pair_id}_snpEff_summary.html
    """
}
