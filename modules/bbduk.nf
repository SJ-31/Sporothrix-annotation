process BBDUK {
    memory { 2.GB * task.attempt }
    publishDir "$outdir/", pattern: "*.txt", mode: 'copy'
    publishDir "$readdir/$name/", pattern: "*B-*fastq.gz", mode: 'copy'
    publishDir "$outdir/", pattern: "*-flagged*"

    input:
    tuple val(name), path(reads) // You NEED the comma
    path(ref)
    val(args)
    val(outdir)
    val(readdir)
    //
    output:
    tuple(val(name), path("B-*.fastq.gz"), emit: reads)
    path("*-flagged*"), emit: flagged
    path("*.txt"), emit: log

    exec:
    def forw = reads[0].baseName.replaceAll(/T-/, 'B-')
    def rev = reads[1].baseName.replaceAll(/T-/, 'B-')
    //
    script:
    """
    bbduk.sh in1=${reads[0]} in2=${reads[1]} \
    out1=${forw}.fastq.gz \
    out2=${rev}.fastq.gz \
    ref=${ref} stats=${name}_bbduk.txt \
    outm=${name}-flagged.fastq.gz  \
    ${args} \
    """
    // bbduk can run out of memory if you don't specify it
}
