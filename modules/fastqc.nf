

process FASTQC {
    publishDir "$", mode: 'symlink'

    input:
    tuple val(id), path(reads)
    val(outdir)
    //
    output:

    //
    script:
    """

    """
    //
}
