/*
 * Module imports
 */
include { MAKER_F; MAKER_GET_FASTA; MAKER_MERGE_FASTA;
    MAKER_COMBINE } from "../modules/maker"
include { MAKER_R as MAKER_R2 } from "../modules/maker"
include { MAKER_R as MAKER_R3 } from "../modules/maker"
include { SNAP } from "../modules/snap"
include { SNAP as SNAP2 } from "../modules/snap"
include { GENEMARKS_ES } from "../modules/genemarksES"
include { AUGUSTUS; AUGUSTUS_MAKER } from "../modules/augustus"
include { BUSCO } from "../modules/busco"
include { FIND_BUSCO } from "../modules/find_buscos"

/*
 * Workflow
 */


workflow train_genemarks {
    take:
    scaffolds

    main:
    genemark_ch = GENEMARKS_ES(scaffolds, params.genemarks)
}

def makerNext = branchCriteria {
            name: !(it =~ /\./)
            fasta: it =~ /fasta/
            gff: it =~ /gff/
            snap: it =~ /hmm/
}

workflow annotation {
    take:
    chromosomes
    scaffolds

    main:
    if ( params.annotate_what == 'scaffolds')
        scaffolds.set { annotating }
    else if ( params.annotate_what == 'chromosome' )
        chromosomes.set { annotating }
    MAKER_F(annotating, params.makerR1, params.outdirAnnotate)
        .set { makerR1 }

    // Train snap
    SNAP(makerR1.gff, '1', params.outdirAnnotate)
        .set { snap1_ch }
    snap1_ch.hmm.join(makerR1.gff).join(annotating).flatten().branch{
            name: !(it =~ /\./)
            fasta: it =~ /fasta/
            gff: it =~ /gff/
            snap: it =~ /hmm/
    }.set { r2_ch }

    // Second round
    MAKER_R2(r2_ch.name, r2_ch.gff, r2_ch.fasta, r2_ch.snap, '2',
    params.makerR2, params.outdirAnnotate)
        .set { makerR2 }
    //  Second SNAP training
    SNAP2(makerR2.gff, '2', params.outdirAnnotate)
        .set { snap2_ch }
    snap2_ch.hmm.join(makerR2.gff).join(annotating).flatten().branch(makerNext)
        .set { r3_ch }

    // Third round
    MAKER_R3(r3_ch.name, r3_ch.gff, r3_ch.fasta, r3_ch.snap, '3',
    params.makerR3, params.outdirAnnotate)
        .set { maker_final }
    MAKER_COMBINE(maker_final.gff.map { it -> [ it[0].replaceAll(/_.*/,''),
    it[1] ]}.groupTuple(), params.outdirAnnotate)

    // Collect results
    MAKER_GET_FASTA(maker_final.makerout, params.outdirAnnotate)
        .set { fastas }
    fastas.protein.map {it -> [ it[0].replaceAll(/_.*/, ''), it[1] ] }
    .transpose()
    .groupTuple().map { it -> [ it[0], 'protein', it[1] ]}
        .set { all_prot }
    fastas.transcripts.map {it -> [ it[0].replaceAll(/_.*/, ''), it[1] ] }
    .transpose()
    .groupTuple().map { it -> [ it[0], 'transcripts', it[1] ]}
        .set { all_transcripts }
    MAKER_MERGE_FASTA(all_prot.mix(all_transcripts), params.outdirAnnotate)
        .set { merged }

    // Verify transcripts
    BUSCO(merged.filter { it[1] =~/transcripts/ }, 'transcriptome', params.outdirAnnotate)
}

workflow get_buscos {
    take:
    busco_dir

    main:
    FIND_BUSCO(busco_dir)
    // Unsure now which busco genes to get
}
