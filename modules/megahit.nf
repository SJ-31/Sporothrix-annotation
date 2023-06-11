process MEGAHIT {
    conda '/home/sc31/Bio_SDD/miniconda3/envs/megahit'

    input:
    tuple val(run), path(reads)
    //
    output:
    path(run)
    //
    script:
    """
    megahit \
    -o $run \
    -1 ${reads[0]} -2 ${reads[1]} \
    """
    //
}

process EXTRACT_MH {
    publishDir "$outdir", mode: 'copy'

    input:
    path(run)
    val(outdir)
    //
    output:
    tuple val("${run}"), path("${run}_megahit.fasta")
    //
    script:
    """
    cp $run/final.contigs.fa ./${run}_megahit.fasta
    """
    //
}
