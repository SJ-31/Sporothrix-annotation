process EXTRACTCHR  {
    publishDir "$outdir/$name", mode: 'copy'
    debug true

    input:
    tuple val(name), path(scaffolds)
    path(ref)
    val(outdir)
    //
    output:
    path("*.fasta")
    //
    shell:
    n_name = scaffolds.baseName.replaceAll(/_.*/, '')
    '''
    rename=!{scaffolds.baseName}
    minimap2 -a !{ref} !{scaffolds} > ${rename}.sam
    samtools sort ${rename}.sam -o ${rename}.bam
    samtools index ${rename}.bam
    samtools view ${rename}.bam | \
    cut -f 3 | sort | uniq > alignments.txt
    while read aln
        do
            case $aln in
            '*')
                samtools view -h ${rename}.bam "${aln}" | \
                samtools fasta > !{n_name}_unaligned.fasta;;
            *)
                samtools view -h ${rename}.bam "${aln}" | \
                samtools fasta > !{n_name}_${aln}.fasta;;
            esac
    done < alignments.txt
    find -type f -empty -delete
    '''
    //
}
