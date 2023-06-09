process FILTERVARIANTS {

    input:
    tuple val(pair_id), path(snps), path(snp_index)
    tuple val(pair_id), path(indels), path(indel_index)
    tuple val(name), path(ref), path(other_ref)
    val(round)

    output:
    tuple (val(pair_id),
    path("${pair_id}_bqsr_snps.vcf"),
    path("${pair_id}_bqsr_indels.vcf"))

    script:
    """
    gatk SelectVariants \
        --exclude-filtered \
        -V $snps \
        -O ${pair_id}_bqsr_snps.vcf
    gatk SelectVariants \
        --exclude-filtered \
        -V $indels \
        -O ${pair_id}_bqsr_indels.vcf
    """
}
