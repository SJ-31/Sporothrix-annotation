process MAFFT {
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
    mafft --preservecase --auto \
    ${combined_fastas} > ${name}_aligned.fasta
    """
}