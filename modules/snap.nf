process SNAP {
    publishDir "$outdir/snap_logs/${sample}_R${round}/", mode: 'copy', pattern: "*.log"

    input:
    tuple val(name), path(gff)
    val(round)
    val(outdir)
    //
    output:
    tuple val(name), path("1-*snap_hmm"), emit: hmm
    path("*.log")
    //
    script:
    sample = name.replaceAll(/_.*/, '')
    """
    maker2zff $gff
    fathom -categorize 100 *.ann *.dna
    fathom genome.ann genome.dna -gene-stats > ${name}_SNAP-gene-stats.log 2>&1
    fathom genome.ann genome.dna -validate > ${name}_SNAP-validate.log 2>&1
    fathom genome.ann genome.dna -categorize 100 > ${name}_SNAP-categorize.log 2>&1
    fathom -export 100 -plus uni.*
    forge export.ann export.dna
    hmm-assembler.pl ${gff.baseName[2..-1]} . > 1-${gff.baseName[2..-1]}.snap_hmm
    """
    //
}

// Alternatively, you could export only 'confident' gene models using the -x <AED threshold> and -l <length threshold> flags in maker2zff


