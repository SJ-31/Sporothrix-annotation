process SNAP {
    publishDir "", mode: 'symlink'

    input:
    tuple val(name), path(gff)
    //
    output:
    tuple val(name), path(snap_hmm)
    //
    script:
    """
    maker2zff $gff
    fathom -categorize 100 *.ann *.dna
    fathom -export 1000 -plus uni.*
    forge export.ann export.dna
    hmm-assembler.pl
    """
    //
}
