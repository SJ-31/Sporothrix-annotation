process SNAP {
    publishDir "$outdir/snap_logs/${sample}_R${round}/", mode: 'copy', pattern: "*.log"
    publishDir "$outdir/snap_models/${sample}_R${round}/", mode: 'copy', pattern: "*.snap_hmm"

    input:
    tuple val(name), path(gff)
    val(round)
    val(outdir)
    //
    output:
    tuple val(name), path("*snap_hmm"), emit: hmm
    path("*.log")
    //
    script:
    sample = name.replaceAll(/_.*/, '')
    """
    maker2zff $gff
    fathom genome.ann genome.dna -gene-stats > ${name}_SNAP-gene-stats.log 2>&1
    fathom genome.ann genome.dna -validate > ${name}_SNAP-validate.log 2>&1
    fathom -categorize 100 *.ann *.dna
    fathom -export 100 -plus uni.*
    forge export.ann export.dna
    hmm-assembler.pl ${gff.baseName} . > ${gff.baseName}.snap_hmm
    """
    //
}
