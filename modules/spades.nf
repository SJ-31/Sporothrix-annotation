process SPADES {
    input:
    tuple val(run), path(reads)
    //
    output:
    tuple val("spades"), val(run), path(run)
    //
    script:
    """
    spades.py \
    -o $run \
    --pe-1 1 ${reads[0]} --pe-2 1 ${reads[1]} \
    """
    //
}

process EXTRACT_SPADES {
    publishDir "$outdir", mode: 'symlink'

    input:
    tuple val(s), val(run), path(assembly)
    val(outdir)
    //
    output:
    tuple val("spades"), path("${run}_scaffolds.fasta")
    //
    script:
    """
    cp $run/scaffolds.fasta ./${run}_scaffolds.fasta
    """
    //
}
