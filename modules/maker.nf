process MAKER_F {
    publishDir "$outdir/makerAllOutput/R1/$sample", pattern: "*.output", mode: 'copy'
    publishDir "$outdir/makerAllOutput/R1/$sample", pattern: "${name}.log", mode: 'copy'
    publishDir "$outdir/${sample}/R1/", mode: 'copy', pattern: "*.gff"
    tag "Annotating $sample, file: $files"

    input:
    tuple val(name), path(files)
    val(args)
    val(outdir)

    output:
    tuple val(name), path("*.output"), emit: makerout
    tuple val(name), path("*.all.gff"), emit: gff
    path("${name}.log")
    //
    script:
    train_on = "opts_genome=${files}"
    sample = name.replaceAll(/_.*/, '')
    predictors = "opts_gmhmm=${args[1]}/${sample}_gmhmm.mod"
    """
    maker -CTL
    maker_cli.py ${args[0]} \
    ${train_on} \
    ${predictors}
    maker -fix_nucleotides 1> ${name}.log
    gff3_merge -d *output*/*datastore*.log
    """
    //
}

process MAKER_R {
    publishDir "$outdir/makerAllOutput/R${round}/$sample", mode: 'copy', pattern: "*all.output"
    publishDir "$outdir/makerAllOutput/R${round}/$sample", mode: 'copy', pattern: "${name}.log"
    publishDir "$outdir/${sample}/R${round}/", mode: 'copy', pattern: "*.gff"
    tag "Round $round: Annotating $sample, gff = $all_gff, snap = $snap, sequence = $fasta"

    input:
    val(name)
    path(all_gff)
    path(fasta)
    path(snap)
    val(round)
    val(args)
    val(outdir)

    output:
    tuple val(name), path("*.output"), emit: makerout
    path("${name}.log")
    tuple val(name), path(all_gff), emit: gff
    //
    script:
    sample = name.replaceAll(/_.*/, '')
    prev = all_gff.baseName.replaceAll(/all.gff/, 'prev.gff')
    """
    maker -CTL
    maker_cli.py ${args} \
    opts_genome=${fasta} \
    opts_maker_gff=${all_gff} \
    opts_snaphmm=${snap}
    maker -fix_nucleotides 1> ${name}.log
    mv $all_gff $prev
    gff3_merge -d *output*/*datastore*.log
    """
    //
}

// Not needed
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
    path("*repeats.gff"), emit: repeats
    path("*protein2genome.gff"), emit: protein
    path("*est2genome.gff"), emit: ests
    tuple val(formatted), path("*{.fasta,.txt}"), emit: fastas

    exec:
    formatted = maker_out.baseName[2..-1].replaceAll(/.maker/, '')
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
