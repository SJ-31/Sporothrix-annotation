process MUSCLE {
    publishDir "$outdir", mode: 'copy'
    input:
    path(combined_fastas)
    val(outdir)
    //
    output:
    path("${group}_*.fasta")
    //
    shell:
    '''
    for fasta in *combined.fasta
        do
            name=$(echo $fasta | sed 's/_.*//')
        muscle -align $fasta -output ${name}_aligned.fasta
        done
    '''
    // Doing it this way is very slow, but necessary if you are running this on a laptop and don't want to run out of memory. Nextflow's queueSize config didn't work for me
}