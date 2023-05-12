process SOAPDENOVO {
    input:
    tuple val(run), path(reads)
    //
    output:

    //
    script:
    """
    cat -e \
    'max_red_len=150' \
    '[LIB]' \
    'avg_ins=144' \
    'reverse_seq=0'
    'q1=${reads[0]}' \
    'q2=${reads[1]}' \
    > soap.config

    """
    //
}
