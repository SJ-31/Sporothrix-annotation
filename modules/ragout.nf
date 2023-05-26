process RAGOUT {
    conda "/home/sc31/Bio_SDD/miniconda3/envs/ragout"
    publishDir "$outdir/ragout"

    input:
    tuple val(name), path(contigs)
    path(reference)
    val(outdir)
    //
    output:
    path("${name}_rag")
    //
    script:
    """
    chef_ragout.py $contigs $reference
    ragout recipe.rcp \
    -o ${name}_rag
    """
    //
}

process EXTRACT_RAGOUT {
    publishDir "$outdir/ragout", mode: 'copy'

    input:
    tuple val(name), path(scaffolds)
    val(outdir)
    //
    output:
    path("${name}.fasta")
    //
    script:
    """
    cp scaffolds/${name}_scaffolds.fasta .
    """
    //
}
