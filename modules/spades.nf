process SPADES {

    input:
    tuple val(run), path(reads)
    //
    output:

    //
    script:
    """
    spades.py \
    -o \
    --pe-1 ${reads[0]} --pe-2 ${reads[1]} \

    """
    //
}
