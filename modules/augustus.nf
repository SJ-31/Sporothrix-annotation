process AUGUSTUS {
    input:
    tuple val(name), path(files)
    val(round)
    //
    output:
    tuple val(name), path("2-$name")
    //
    script:
    """
    touch 2-$name
    """
    //
}
