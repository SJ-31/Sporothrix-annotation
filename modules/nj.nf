process NJ {
    publishDir "$outdir/scores", mode: 'copy', pattern: "*score.txt"
    publishDir "$outdir/trees", mode: 'copy', pattern: "*NJ.nwk"

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