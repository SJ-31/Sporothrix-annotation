process FILTER_MT {
    publishDir "$outdir", mode: 'copy', pattern: "*mtDNA.fasta"

    input:
    tuple val(name), path(scaffolds)
    path(mtdna)
    val(outdir)
    //
    output:
    tuple val(name), path("${scaffolds}"), emit: filtered
    path("*mtDNA.fasta")
    //
    script:
    """
    mv ${scaffolds} unfiltered.fasta
    minimap2 -k 28 \
    unfiltered.fasta ${mtdna} > alignment.paf
    makeblastdb -in unfiltered.fasta -parse_seqids -blastdb_version 5 \
    -out blastdb -dbtype nucl
    blastn -query $mtdna -db blastdb -outfmt 6 > ${name}_blast.txt
    cut -f 6 alignment.paf > loc.txt
    seqkit grep -v -f loc.txt unfiltered.fasta > ${scaffolds.baseName}.fasta
    seqkit grep -f loc.txt unfiltered.fasta > ${name}_mtDNA.fasta
    """
    // Note the reversed order here: we use the scaffolds as the reference to find mtDNA
}
