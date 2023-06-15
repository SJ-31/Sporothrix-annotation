process GFF_BUSCO_MAP {
    publishDir "$outdir", mode: 'copy', pattern: "*tsv"

    input:
    tuple val(name), path(liftover), path(fasta)
    val(busco_gff_ref)
    val(outdir)
    //
    output:
    tuple val(name), path("${name}_busco_gff.tsv"), path(fasta)
    //
    script:
    """
    liftoff_transfer.sh $name $liftover $busco_gff_ref ${name}_liftoff_busco.tsv
    """
    //
}
// You need to use the busco_gff file generated on the original reference sequence to translate the new locations of the busco genes on the gff lifted over from the sample