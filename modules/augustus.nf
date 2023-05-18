process AUGUSTUS {
    input:
    tuple val(name), path(file)
    val(round)
    //
    output:
    tuple val(name), path("2-${round}_${name}")
    //
    exec:
    def species = "${round}_${name}"
    script:
    """
    touch 2-$species
    """
    //
}

process AUGUSTUS_MAKER {
    publishDir "$outdir/$name/$round", mode: 'copy', pattern: "*{.out,log}"
    input:
    tuple val(name), path(file)
    val(round)
    val(outdir)
    val(augustus_config)
    //
    output:
    tuple val(name), path("2-${round}_${name}")
    path("*.out")
    //
    exec:
    def assembly = file[0]
    def gff = file[1]
    def initialGB = "${name}.gb"
    def evaluate = "${name}.gb.evaluate"
    def species = "${round}_${name}"

    script:
    """
    augustus_maker.sh $assembly $gff $name $initialGB $evaluate $species
    """
    //
}

process AUGUSTUS_BUSCO {
    input:

    //
    output:

    //
    script:
    """

    """
    //
}
