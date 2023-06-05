process FIND_BUSCO {
    publishDir "$outdir/$name", mode: "copy"

    input:
    tuple val(name), path(busco_dir)
    val(params.wanted_buscos)
    val(outdir)
    //
    output:
    tuple val(name), path("*found*")
    //
    script:
    """
    get_busco.py $busco_dir \
    ${name}-selected_buscos.fasta \
    $params.wanted_buscos
    """
    //
}
