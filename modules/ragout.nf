process RAGOUT {
    conda "/home/sc31/Bio_SDD/miniconda3/envs/ragout"

    input:
    tuple val(name), path(contigs)
    path(reference)
    val(outdir)
    //
    output:
    path("${n_name}_rag")
    exec:
    n_name = name.replaceAll(/_.*/, '').replaceAll(/-.*/)
    //
    script:
    """
    chef_ragout.py $contigs $reference
    ragout recipe.rcp \
    -o ${n_name}_rag
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
    name = scaffolds.baseName.replaceAll(/_.*/, '').replaceAll(/-.*/)

    script:
    """
    cp $scaffolds/*_scaffolds.fasta .
    """
    //
}
