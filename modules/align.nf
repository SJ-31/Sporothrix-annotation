process ALIGN {
    tag "Aligning with $ref"

    input:
    tuple val(pair_id), path(reads)
    tuple val(name), path(ref), path(other_ref)

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
