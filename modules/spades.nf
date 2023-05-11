process SPADES {
    publishDir "$outdir", mode: 'symlink', pattern: "$run"

    input:
    tuple val(run), path(reads)
    val(outdir)
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
    input:
    path(run)
    //
    output:

    //
    script:
    """

    """
    //
}
