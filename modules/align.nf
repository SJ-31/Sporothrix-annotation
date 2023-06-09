process ALIGN {
    tag "Aligning with $ref"
    publishDir "${outdir}/aligned_reads/$pair_id", mode:'copy'

    input:
    tuple val(pair_id), path(reads)
    tuple val(name), path(ref), path(other_ref)
    val(outdir)

    output:
    tuple val(pair_id), path("${pair_id}_aligned_reads.sam") \

    script:
    """
    bwa mem \
	-K 100000000 \
	-v 3 \
	-t ${task.cpus} \
	-Y \
	$ref \
	${reads[0]} \
	${reads[1]} \
	> ${pair_id}_aligned_reads.sam
    """
}
