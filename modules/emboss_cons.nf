process CONS {
    input:
    path(alignment)
    val(outdir)
    //
    output:
    val("*cons.fasta")
    //
    shell:
    '''
    header=$(sed 's/>S[[:digit:]]\{2\}|/>cons|/g')
    name=$(echo !{alignment} | sed 's/_.*//g')
    cons -sequence !{alignment} \
        -outseq ${name}_cons.fasta \
        -name $header
    '''
    //
}