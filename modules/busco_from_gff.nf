process BUSCO_FROM_GFF {

    input:
    tuple val(name), path(busco_gff_tsv), path(gff), path(fasta)
    //
    output:
    tuple val(name), path("*.fasta")
    //
    script:
    """
    busco_from_gff.sh -s $name -b $busco_gff_tsv -g $gff -f $fasta
    """
    //
}
