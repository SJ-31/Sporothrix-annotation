process NJ {
    publishDir "$outdir", mode: 'copy', pattern: "*score.txt"

    input:
    path(alignment)
    val(reference)
    val(outdir)
    //
    output:
    tuple path("*score.txt"), path("*_NJ.nwk")
    //
    script:
    """
    Rscript $bin/neighbor_joining.r $reference $alignment
    """
    //
}