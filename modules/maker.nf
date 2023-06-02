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
    publishDir "$outdir/makerAllOutput/R${round}/$sample", mode: 'copy', pattern: "*maker.output"
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

process COMBINE {
    publishDir "$outdir/${sample}/", mode: 'copy', pattern: "*final.gff"

    input:
    tuple val(sample), path(gffs)
    val(outdir)
    //
    output:
    path("${sample}_final.gff")
    //
    script:
    """
    cat $gffs > ${sample}_final.gff
    """
    //
}

process GET_FASTA {
    publishDir "$outdir/${sample}/protein_fastas/", mode: 'copy', pattern: "*protein*"
    publishDir "$outdir/${sample}/transcript_fastas/", mode: 'copy', pattern: "*transcripts*"

    input:
    tuple val(name), path(maker_out)
    val(outdir)
    //
    output:
    tuple val(name), path("*protein*"), emit: protein
    tuple val(name), path("*transcripts*"), emit: transcripts

    script:
    sample = name.replaceAll(/_.*/, '')
    """
    fasta_merge -d *output*/*datastore*.log
    """
    //
}

process MERGE_FASTA {
    publishDir "$outdir/${sample}/", mode: 'copy', pattern: "*.fasta"

    input:
    tuple val(sample), path(files)
    val(outdir)
    //
    output:
    path("${sample}_merged.fasta")

    script:
    sample = name.replaceAll(/_.*/, '')
    """
    cat $files > ${sample}_merged.fasta
    """
    //
}
