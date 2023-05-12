process MAKER {
    input:

    //
    output:

    //
    script:
    """
    maker -CTL
    maker_cli.py \
    est2genome=1 \
    protein2genome=2
    maker
    """
    //
}
