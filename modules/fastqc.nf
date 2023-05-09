

process FASTQC {
    publishDir "$outdir/$run", pattern: "*.html", mode: 'copy'

    input:
    tuple val(run), path(reads)
    val(outdir)
    //
    output:
    path("*.html"), emit: html
    path("*.zip"), emit: zip
    //
    script:
    """
    fastqc $reads
    """
    //
}
