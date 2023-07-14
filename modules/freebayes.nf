process FREEBAYES {
    publishDir "$outdir/variants/$pair_id", mode: 'copy'

    input:
    path(alignment)
    path(reference)
    val(outdir)
    //
    output:
    tuple (val(pair_id), path("${pair_id}_merged.vcf.gz"))
    //
    script:
    """
    samtools view -S -b $alignment > alignment.bam
    freebayes -f $reference alignment.bam > ${pair_id}_merged.vcf.gz
    """
    //
}
