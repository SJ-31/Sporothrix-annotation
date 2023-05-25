process SSPACE {
    publishDir "$outdir/sspace/$sample", mode: 'copy'

    input:
    tuple val(sample), path(contigs), path(reads)
    val(outdir)
    val(args)
    //
    output:
    tuple val(sample), path("*final*")
    //
    script:
    """
    echo -n "lib1" > lib.txt
    echo -n " ${reads[0]} ${reads[1]} " >> lib.txt
    echo -n "151 0.5 FR" >> lib.txt
    SSPACE_Basic_v2.0.pl \
    $args \
    -l lib.txt \
    -s $contigs
    """
    //
}
