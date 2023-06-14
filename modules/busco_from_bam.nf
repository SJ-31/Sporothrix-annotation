process BUSCO_FROM_BAM {
    debug true

    input:
    tuple val(name), path(sample)
    val(reference)
    val(busco_ref)
    //
    output:
    path("*at*fasta")
    //
    script:
    """
    minimap2 -a $reference $sample > aligned.sam
    samtools sort aligned.sam -o aligned.bam
    samtools index aligned.bam
    bfb_busco_gff.sh $name $busco_ref aligned.bam
    """
    // fields[0] = Gene busco ID
    // fields[1] = Gene sequence location
    // fields[2], fields[3] = Gene start, stop
    // fields[4] = Gene strandedness
    // fields[5] = Gene description
}