process MUSCLE {
    publishDir "$outdir", mode: 'copy'
    input:
    path(combined_fastas)
    val(outdir)
    //
    output:
    path("*aligned.fasta")
    //
    script:
    name = combined_fastas.baseName.replaceAll(/_.*/, '')
    """
    muscle -align ${combined_fastas} -output ${name}_aligned.fasta
    """
}