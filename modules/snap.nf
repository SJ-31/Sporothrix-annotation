process SNAP {
    publishDir "$params.outdirAnnotate/$name/$round/", mode: 'symlink', pattern: "*.log"

    input:
    tuple val(name), path(gff)
    val(round)
    //
    output:
    tuple val(name), path("1-*snap_hmm"), emit: hmm
    path("*.log")
    //
    script:
    """
    maker2zff $gff
    fathom -categorize 100 *.ann *.dna
    fathom genome.ann genome.dna > SNAP-gene-stats.log 2>&1
    fathom genome.ann genome.dna > SNAP-validate.log 2>&1
    fathom -export 1000 -plus uni.*
    forge export.ann export.dna
    hmm-assembler.pl ${gff.baseName[2..-1]} . > 1-${gff.baseName[2..-1]}.snap_hmm
    """
    //
}

// Alternatively, you could export only 'confident' gene models using the -x <AED threshold> and -l <length threshold> flags in maker2zff


