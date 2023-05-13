process MULTIQC {
    publishDir "$outdir/$round", pattern: "*html", mode: 'copy'
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
