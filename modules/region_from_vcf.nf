process VCF_GET_REGION {
    publishDir "$outdir/$name", mode: 'copy', pattern: "${name}*.vcf"
    publishDir "$outdir", mode: 'copy', pattern: "*num_vars.tsv"

    input:
    tuple val(name), path(vcf)
    val(outdir)
    val(gene_locs)
    //
    output:
    path("${name}_*.vcf"), emit: vars
    path("*num_vars.tsv*"), emit: nums
    //
    script:
    """
    bcftools index $vcf
    all_vars.sh $name $gene_locs $vcf
    """
    //
}