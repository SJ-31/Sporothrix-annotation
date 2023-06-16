process KALLISTO {
    input:
    tuple val(sample), path(combined), path(reads)
    //
    output:

    //
    script:
    """
    kallisto index -i index 
    """
    //
}