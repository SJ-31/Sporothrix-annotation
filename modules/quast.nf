process QUAST {
    publishDir "$outdir/$assembler", mode: 'copy'

    input:
    tuple val(assembler), path(assemblies)
    val(outdir)
    path(ref_genome)
    path(ref_genes)
    //
    output:

    //
    script:
    """
    quast.py $assemblies \
    --fungal \
    -r $ref_genome \
    -g $ref_genes
    """
    //
}
