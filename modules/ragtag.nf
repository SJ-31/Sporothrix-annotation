process RAGTAG {
    conda "/home/sc31/Bio_SDD/miniconda3/envs/ragtag"
    publishDir "$outdir/ragtag", mode: 'copy'

    input:
    tuple val(name), path(contigs)
    path(reference)
    val(outdir)
    //
    output:
    path("${name}-ragtag_scaffolds.fasta")
    //
    script:
    """
    ragtag.py scaffold $reference $contigs $reference
    cp ragtag_output/ragtag.scaffold.fasta ./${name}-ragtag_scaffolds.fasta
    """
    //
}
