process GENEMARKS_ES {
    input:
    tuple val(name), path(file)
    //
    output:
    tuple val(name), path("1-gmhmm.mod")
    //
    script:
    """
    gmes_petap.pl -ES -sequence ${file} -fungus
    mv gmhmm.mod 1-gmhmm.mod
    """
    //
}
