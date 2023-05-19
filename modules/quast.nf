process QUAST {
    publishDir "$outdir/$assembler", mode: 'copy'

    input:
    path(assemblies)
    val(outdir)
    path(ref_genome)
    path(ref_genes)
    //
    output:
    path("*")
    //
    script:
    """
    quast.py $assemblies \
    --fungus \
    -r $ref_genome \
    -g $ref_genes
    """
    //
}
