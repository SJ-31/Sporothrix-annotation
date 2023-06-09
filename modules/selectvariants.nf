process SELECTVARIANTS {
    input:
    tuple val(pair_id), file(raw_variants)
    tuple val(name), path(ref), path(other_ref)
    val(round)

    output:
    tuple val(pair_id), path("${pair_id}_raw_snps_${round}.vcf"), emit: raw_snps_ch
    tuple val(pair_id), path("${pair_id}_raw_indels_${round}.vcf"), emit: raw_indels_ch

    script:
    """
    gatk SelectVariants \
	-R $ref \
	-V $raw_variants \
	-select-type SNP \
	-O ${pair_id}_raw_snps_${round}.vcf
    gatk SelectVariants \
        -R $ref \
        -V $raw_variants \
        -select-type INDEL \
        -O ${pair_id}_raw_indels_${round}.vcf
    """
}
