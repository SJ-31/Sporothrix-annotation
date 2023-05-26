process RAGTAG {
    conda "/home/sc31/Bio_SDD/miniconda3/envs/ragtag"
    publishDir "$outdir/ragtag"

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
    path("*.fasta")
    //
    script:
    def name = scaffolds.baseName.replaceAll(/_.*/, '')
    """
    cp $scaffolds/*_scaffolds.fasta .
    """
    //
}
