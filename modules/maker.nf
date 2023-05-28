process MAKER_F {
    input:
    tuple val(sample), path(files)
    val(args)

    output:
    tuple val(sample), path("*.output")
    //
    script:
    train_on = "opts_genome=${files}"
    predictors = "opts_gmhmm=${args[1]}/1-${sample}_gmhmm.mod"
    """
    maker -CTL
    maker_cli.py ${args[0]} \
    ${train_on} \
    ${predictors}
    maker -fix_nucleotides
    """
    //
}

process MAKER_R {
    input:
    tuple val(sample), path(files)
    val(args)

    output:
    tuple val(sample), path("*.output")
    //
    script:
    train_on = "opts_genome=${files}"
    predictors = """
        opts_snaphmm=${files[1]} \
        opts_est_gff=${files[2]} \
        opts_protein_gff=${files[3]} \
        opts_rm_gff=${files[4]}
        """
        // genemarks = "opts_gmhmm=$current[1]" // Maybe won't need genemarks in the second round since it hasn't been trained on anything else
    """
    maker -CTL
    maker_cli.py ${args[0]} \
    ${train_on} \
    ${predictors}
    maker -fix_nucleotides
    """
    //
}

process GET_GFF {
    publishDir "$outdir/${name}/R${round}/", mode: 'copy', pattern: "*.all.gff"
    publishDir "$outdir/${name}/R${round}/", mode: 'copy', pattern: "*.fasta"
    publishDir "$outdir/${name}/R${round}/", mode: 'copy', pattern: "*.txt"
    debug true
    input:
    tuple val(name), path(maker_out)
    val(round)
    val(outdir)
    //
    output:
    tuple val(formatted), path("0-*.all.gff"), emit: all
    tuple val(formatted), path("*{est2genome,repeats,protein2genome}.gff"), emit: evidence
    tuple val(formatted), path("*{.fasta,.txt}"), emit: fastas

    exec:
    formatted = maker_out.baseName[2..-1].replaceAll(/.maker/, '')

    shell:
    '''
    name=$(echo *output* | sed -e 's/0-//' -e 's/.maker.output//' )
    cp -r *output/* .
    gff3_merge -d 0-${name}_master_datastore_index.log
    fasta_merge -d 0-${name}_master_datastore_index.log
    gff3_merge -n -s -d 0-${name}_master_datastore_index.log \
    > ${name}.all.maker.noseq.gff
    awk '{ if ($2 == "est2genome") print $0 }' \
    ${name}.all.maker.noseq.gff > 2-${name}.all.maker.est2genome.gff
    awk '{ if ($2 == "protein2genome") print $0 }' \
    ${name}.all.maker.noseq.gff > 3-${name}.all.maker.protein2genome.gff
    awk '{ if ($2 ~ "repeat") print $0 }' \
    ${name}.all.maker.noseq.gff > 4-${name}.all.maker.repeats.gff
    if [ ! -f ./*fasta ];
        then
            touch 0-${name}.all.NOFASTA.txt
    fi
    '''
    //
}
