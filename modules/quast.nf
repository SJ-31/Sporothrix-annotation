process QUAST {
    publishDir "$outdir/", mode: 'copy'

    input:
    path(assemblies)
    val(outdir)
    path(ref_genome)
    path(ref_genes)
    val(args)
    //
    output:
    path("*")
    //
    script:
    """
    quast.py $assemblies \
    --fungus \
    -r $ref_genome \
    -g $ref_genes \
    $args
    """
    //
}
