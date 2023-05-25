process REPEATMODELER {
    input:
    tuple val(species), path(genome)
    //
    output:
    tuple val(species), path("*")
    //
    script:
    """
    BuildDatabase -name ${species}_db $genome
    RepeatModeler -database ${species}_db -LTRStruct
    """
    //
}

process GET_REPEAT {
    publishDir "$outdir/$species", mode: 'copy'

    input:
    tuple val(species), path(rm)
    val(outdir)
    //
    output:
    path("*")
    //
    script:
    """
    cp -r RM*/* .
    """
    //
}
