process MAKER_F {
    publishDir "$outdir/makerAllOutput/R1/$sample", pattern: "*.output", mode: 'copy'
    publishDir "$outdir/${sample}/R1/", mode: 'copy', pattern: "*.gff"
    tag "Annotating $sample, file: $files"

    input:
    tuple val(sample), path(files)
    val(args)
    val(outdir)

    output:
    tuple val(sample), path("*.output"), emit: makerout
    tuple val(sample), path("*.all.gff"), emit: all
    //
    script:
    train_on = "opts_genome=${files}"
    predictors = "opts_gmhmm=${args[1]}/${sample}_gmhmm.mod"
    """
    maker -CTL
    maker_cli.py ${args[0]} \
    ${train_on} \
    ${predictors}
    maker -fix_nucleotides \
    gff3_merge -d ${sample}_master_datastore_index.log
    """
    //
}

process MAKER_R {
    publishDir "$outdir/$round/$name", mode: 'copy', pattern: "*all.output"
    publishDir "$outdir/${name}/R${round}/", mode: 'copy', pattern: "*.gff"
    tag "Round $round: Annotating $name in $prevout"

    input:
    tuple val(name), path(prevout)
    val(snap)
    val(round)
    val(args)
    val(outdir)

    output:
    tuple val(name), path("*.output"), emit: makerout
    tuple val(name), path("*all.gff"), emit: gff
    //
    script:
    """
    maker -CTL
    maker_cli.py ${args} \
    maker -fix_nucleotides \
    -base ${name}_R${round}
    gff3_merge_d ${name}_master_datastore_index.log
    """
    //
}

// bug: Doing this way doesn't seem to be working
// todo: Fix this up!!!
process MAKER_B {
    tag "Annotating $chromosome,$fasta,$est2genome,$protein2genome,$repeats,$snap"

    input:
    tuple val(chromosome)
    path(fasta)
    path(est2genome)
    path(protein2genome)
    path(snap)
    path(repeats)
    val(args)

    output:
    tuple val(sample), path("*.output")
    //
    script:
    sample = chromosome.replaceAll(/_.*/, '')
    """
    maker -CTL
    maker_cli.py ${args} \
    opts_genome=$fasta \
    opts_snaphmm=$snap \
    opts_est_gff=$est2genome \
    opts_protein_gff=$protein2genome \
    opts_rm_gff=$repeats
    maker -fix_nucleotides
    """
    //
}

process GET_GFF {
    tag "Extracting $name, round $round"

    publishDir "$outdir/${name}/R${round}/", mode: 'copy', pattern: "*.gff"
    publishDir "$outdir/${name}/R${round}/", mode: 'copy', pattern: "*.fasta"
    publishDir "$outdir/${name}/R${round}/", mode: 'copy', pattern: "*.txt"
    debug true
    input:
    tuple val(name), path(maker_out)
    val(round)
    val(outdir)
    //
    output:
    tuple val(formatted), path("*.all.gff"), emit: all
    tuple val(formatted), path("*{est2genome,repeats,protein2genome}.gff"), emit: evidence
    tuple val(formatted), path("*{.fasta,.txt}"), emit: fastas

    exec:
    formatted = maker_out.baseName[2..-1].replaceAll(/.maker/, '')
    //todo: Something is wrong with the way you are extracting the protein and genome stuff
    shell:
    '''
    name=$(echo *output* | sed -e 's/.maker.output//' )
    cp -r *output/* .
    gff3_merge -d ${name}_master_datastore_index.log
    fasta_merge -d ${name}_master_datastore_index.log
    gff3_merge -n -s -d ${name}_master_datastore_index.log \
    > ${name}.all.maker.noseq.gff
    awk '{ if ($2 == "est2genome") print $0 }' \
    ${name}.all.maker.noseq.gff > ${name}.all.maker.est2genome.gff
    awk '{ if ($2 == "protein2genome") print $0 }' \
    ${name}.all.maker.noseq.gff > ${name}.all.maker.protein2genome.gff
    awk '{ if ($2 ~ "repeat") print $0 }' \
    ${name}.all.maker.noseq.gff > ${name}.all.maker.repeats.gff
    if [ ! -f ./*fasta ];
        then
            touch ${name}.all.NOFASTA.txt
    fi
    '''
    //
}
