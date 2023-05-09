process MULTIQC {
    publishDir "$outdir", pattern: "*html", mode: 'copy'
    input:
    path(combined)
    val(outdir)
    //
    output:
    path("*")
    //
    script:
    """
    multiqc $combined
    """
    //
}
