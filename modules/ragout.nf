process RAGOUT {
    conda "/home/sc31/Bio_SDD/miniconda3/envs/ragout"
    publishDir "$outdir/ragout"

    input:
    tuple val(name), path(contigs)
    path(reference)
    val(outdir)
    //
    output:
    path("${name}_rag") //todo: add path
    //
    script:
    """
    echo ".references = $reference.baseName" > recipe.rcp
    echo ".target = $contigs.baseName\n" >> recipe.rcp
    echo "$contigs = ./$contigs" >> recipe.rcp
    echo "$reference = ./$reference" >> recipe.rcp
    ragout recipe.rcp \
    -o ${name}_rag
    """
    //
}
