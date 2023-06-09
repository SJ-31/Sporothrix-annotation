process ANALYZECOVARIATES {
    publishDir "${outdir}/bqsr/$pair_id", mode:'copy'

    input:
    tuple val(pair_id), path(recal_table), path(post_recal_table)
    val(outdir)

    output:
    tuple val(pair_id), path("${pair_id}_recalibration_plots.pdf")

    script:
    """
    gatk AnalyzeCovariates \
	-before $recal_table \
	-after $post_recal_table \
	-plots ${pair_id}_recalibration_plots.pdf
    """
}
