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
    publishDir "$outdir2", mode: 'copy', pattern: "${species}_RM.fa.classified"

    input:
    tuple val(species), path(rm)
    val(outdir)
    val(outdir2)
    //
    output:
    path("*")
    //
    script:
    """
    cp -r RM*/* .
    mv consensi.fa.classified ${species}_RM.fa.classified
    """
    //
}
