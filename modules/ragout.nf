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
    path(scaffolds)
    val(outdir)
    //
    output:
    tuple val(name), path("*.fasta")
    //
    exec:
    name = scaffolds.baseName.replaceAll(/_.*/, '')

    script:
    """
    cp $scaffolds/*_scaffolds.fasta .
    """
    //
}
