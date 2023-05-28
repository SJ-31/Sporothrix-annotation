process GENEMARKS_ES {
    publishDir "$outdir", pattern: "1-${name}_gmhmm.mod", mode: 'copy'

    input:
    tuple val(name), path(file)
    val(outdir)
    //
    output:
    path("1-${name}_gmhmm.mod")

    //
    script:
    """
    gmes_petap.pl -ES -sequence ${file} -fungus
    mv gmhmm.mod 1-${name}_gmhmm.mod
    """
    //
}
