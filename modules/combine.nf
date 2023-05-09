process COMBINE {
    input:
    path("*")
    //
    output:
    path("combined")
    //
    script:
    """
    mkdir combined
    ls --ignore=combined | xargs -I{} mv {} combined/
    """
    //
}
