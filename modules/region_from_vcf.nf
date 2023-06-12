process VCF_GET_REGION {
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(name), path(vcf)
    val(outdir)
    each name_regions
    //
    output:
    path("${name}_${feature}.vcf")
    //
    script:
    (feature, range) = name_regions.tokenize('|')
    """
    bcftools index $vcf
    bcftools view -r ${range} $vcf > ${name}_${feature}.vcf
    """
    //
}