process SPADES {
    input:
    tuple val(run), path(reads)
    path(args)
    output:
    tuple val("spades"), val(run), path(run)
    //
    script:
    """
    spades.py \
    -o $run \
    $args \
    --pe-1 1 ${reads[0]} --pe-2 1 ${reads[1]}
    """
    //
}

process RNA_SPADES {
    input:
    tuple val(run), path(reads)
    //
    output:
    tuple val("spades"), val(run), path(run)
    //
    script:
    """
    rnaspades.py \
    -o $run \
    --pe-1 1 ${reads[0]} --pe-2 1 ${reads[1]} \
    """
    //
}


process EXTRACT_SPADES {
    publishDir "$outdir", mode: 'copy'

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

process EXTRACT_RNASPADES {
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(s), val(run), path(assembly)
    val(outdir)
    //
    output:
    tuple val("spades"), path("${run}_transcripts.fasta")
    //
    script:
    """
    cp $run/transcripts.fasta ./${run}_transcripts.fasta
    """
    //
}
