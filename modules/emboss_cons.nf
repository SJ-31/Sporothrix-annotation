process CONS {
    publishDir "$outdir", mode: 'copy'

    input:
    path(alignment)
    val(outdir)
    //
    output:
    path("*cons.fasta")

    //
    script:
    """
    emboss_cons.sh $alignment
    """
    //
}