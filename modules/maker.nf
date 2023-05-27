process MAKER {
    input:
    tuple val(name), path(files)
    val(args)
    val(first) // Controls whether or not maker is trained on the genome or the gff file from a previous round
    output:
    tuple val(name), path("*.output")
    //
    script:
    train_on = "opts_genome=${files[0]}"
    if (first) {
    predictors = "opts_gmhmm=${files[1]}"
    }
    else{
        predictors = """
            opts_snaphmm=${files[1]} \
            opts_est_gff=${files[3]} \
            opts_protein_gff=${files[4]} \
            opts_rm_gff=${files[5]}
            """
        // genemarks = "opts_gmhmm=$current[1]" // Maybe won't need genemarks in the second round since it hasn't been trained on anything else
    }
    """
    maker -CTL
    maker_cli.py ${args} \
    ${train_on} \
    ${predictors}
    maker -fix_nucleotides
    """
    //
}

process GET_GFF {
    publishDir "$outdir/$name/$round/", mode: 'copy', pattern: "*.all.gff"
    publishDir "$outdir/$name/$round/", mode: 'copy', pattern: "*.fasta"

    debug true
    input:
    tuple (val(name), (path(maker_out)))
    val(round)
    //
    output:
    tuple val(name), path("0-${name}.all.gff"), emit: all
    tuple val(name), path("*{est2genome,repeats,protein2genome}.gff"), emit: evidence
    tuple val(name), path("*.fasta"), emit: fastas

    shell:
    '''
    cp -r *output/* .
    gff3_merge -d 0-!{name}_master_datastore_index.log
    fasta_merge -d 0-!{name}_master_datastore_index.log
    gff3_merge -n -s -d 0-!{name}_master_datastore_index.log \
    > !{name}.all.maker.noseq.gff
    awk '{ if ($2 == "est2genome") print $0 }' \
    !{name}.all.maker.noseq.gff > 3-!{name}.all.maker.est2genome.gff
    awk '{ if ($2 == "protein2genome") print $0 }' \
    !{name}.all.maker.noseq.gff > 4-!{name}.all.maker.protein2genome.gff
    awk '{ if ($2 ~ "repeat") print $0 }' \
    !{name}.all.maker.noseq.gff > 5-!{name}.all.maker.repeats.gff
    '''
    //
}
