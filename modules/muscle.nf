process MUSCLE {
    publishDir "$outdir", mode: 'copy'
    input:
    tuple val(group), path(samples)
    val(outdir)
    //
    output:
    path("${group}_aligned.fasta")
    //
    script:
    """
    cat $samples > combined.fasta
    muscle -align combined.fasta -output ${group}_aligned.fasta
    """
    //
}