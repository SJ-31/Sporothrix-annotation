process SEQKIT_STATS {
    input:
    tuple val(phase), path(reads)
    //
    output:
    path(formatted)
    //
    script:
    """
    seqkit stats $reads
    -e \
    >  ${phase}_seqkit_stats.txt
    """
    //
}
