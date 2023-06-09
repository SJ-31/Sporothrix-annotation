process HAPLOTYPECALLER {
    input:
    tuple val(pair_id), path(input_bam)
    tuple val(name), path(ref), path(other_ref)
    val(round)

    output:
    tuple val(pair_id), path("${pair_id}_raw_variants_${round}.vcf")

    script:
    """
    gatk HaplotypeCaller \
	-R $ref \
	-I $input_bam \
	-O ${pair_id}_raw_variants_${round}.vcf \
    """
}
