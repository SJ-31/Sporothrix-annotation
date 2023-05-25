process CAP {
    publishDir "$outdir/$name"

    input:
    tuple val(name), path(contigs)
    val(outdir)
    //
    output:
    path("*CAP.fasta")
    //
    script:
    """
    cap3 $contigs
    cat ${contigs}.cap.contigs ${contigs}.cap.singles > ${contigs.baseName}_CAP.fasta
    """
    //
}
