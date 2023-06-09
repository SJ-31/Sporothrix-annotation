process MAKE_REF {
    input:
    path(fasta)
    //
    output:
    tuple val(name), path(fasta), path("*{.fasta.,dict}*")
    //
    script:
    name = fasta.baseName.replaceAll(/-.*/, '').replaceAll(/.*-/, '')
    """
    bwa index $fasta
    samtools faidx $fasta
    java -jar $PICARD_JAR CreateSequenceDictionary \
    R=$fasta \
    O=${fasta.baseName}.dict
    """
    //
}
