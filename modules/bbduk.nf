process BBDUK {
    tag "Cleaning with $reference"
    memory { 2.GB * task.attempt }
    publishDir "$outdir/", pattern: "{*.txt, *-flagged*}", mode: 'copy'

    input:
    tuple val(name), path(reads) // You NEED the comma
    path(ref)
    val(args)
    val(outdir)
    //
    output:
    tuple(val(reference.baseName), path("B-*.fastq.gz"), emit: reads)
    path("*-flagged*")
    path("*.txt")

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
