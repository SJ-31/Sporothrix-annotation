process GFF_BUSCO_MAP {
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(name), path(gff), path(fasta)
    val(busco_ref)
    val(outdir)
    //
    output:
    tuple val(name), path("${name}_busco_gff.tsv"), path(fasta)
    //
    script:
    """

    """
    //
}
// You need to use the busco_gff file generated on the original reference sequence to translate the new locations of the busco genes on the gff lifted over from the sample