process COMBINE {
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(group), path(samples)
    path(outdir)
    //
    output:
    path({"${group}_{combined,dne}.fasta"})
    script:
    """
    cat $samples > ${group}_combined.fasta
    [ -s ${group}_combined.fasta ] || touch ${group}_dne.fasta
    """
}
