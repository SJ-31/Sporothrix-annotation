process COMBINE {
    publishDir "$outdir", mode: 'copy', pattern: "*combined.fasta"
    publishDir "$outdir/MISSING", mode: 'copy', pattern: "*dne*"

    input:
    tuple val(group), path(samples)
    val(outdir)
    //
    output:
    path({"${group}_combined.fasta"}), emit: exists, optional: true
    path("*dne.fasta"), emit: dne, optional: true
    script:
    """
    cat $samples > ${group}_combined.fasta
    if [ -s ${group}_combined.fasta ]; then
        echo "Exists"
    else
        touch ${group}_dne.fasta
        rm ${group}_combined.fasta
    fi
    """
}
