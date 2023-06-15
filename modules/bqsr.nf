process BQSR {

    input:
    tuple (val(pair_id),
	path(input_bam),
	path(filtered_snps),
	path(filtered_indels))
    tuple val(name), path(ref), path(other_ref)

    output:
    tuple val(pair_id), path("${pair_id}_recal_data.table"), path("${pair_id}_post_recal_data.table"), emit: analyze_covariates_in_ch
    tuple val(pair_id), path("${pair_id}_recal.bam"), emit: recalibrated_bams

    script:
    """
    gatk IndexFeatureFile \
        --I $filtered_snps
    gatk IndexFeatureFile \
        --I $filtered_indels

    gatk BaseRecalibrator \
	-R $ref \
	-I $input_bam \
	--known-sites ${pair_id}_bqsr_snps.vcf \
	--known-sites ${pair_id}_bqsr_indels.vcf \
	-O ${pair_id}_recal_data.table

    gatk ApplyBQSR \
        -R $ref \
        -I $input_bam \
        -bqsr ${pair_id}_recal_data.table \
        -O ${pair_id}_recal.bam

    gatk BaseRecalibrator \
        -R $ref \
	-I ${pair_id}_recal.bam \
        --known-sites ${pair_id}_bqsr_snps.vcf \
	--known-sites ${pair_id}_bqsr_indels.vcf \
	-O ${pair_id}_post_recal_data.table
    """
}
