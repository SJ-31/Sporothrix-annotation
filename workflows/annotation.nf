/*
 * Module imports
 */
include { GET_GFF } from "../modules/maker"
include { GET_GFF as GET_GFF2 } from "../modules/maker"
include { GET_GFF as GET_GFF3 } from "../modules/maker"
include { MAKER_F } from "../modules/maker"
include { MAKER_R as MAKER_R2 } from "../modules/maker"
include { MAKER_R as MAKER_R3 } from "../modules/maker"
include { SNAP } from "../modules/snap"
include { GENEMARKS_ES } from "../modules/genemarksES"
include { AUGUSTUS; AUGUSTUS_MAKER } from "../modules/augustus"

/*
 * Workflow
 */

workflow train_genemarks {
    take:
    scaffolds

    main:
    genemark_ch = GENEMARKS_ES(scaffolds, params.genemarks)
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
    // round1_out = GET_GFF(makerR1, '1', params.outdirAnnotate)
    // round1_out.evidence.transpose()
    //     .set { r1_evidence_ch }
    // Train predictors
    SNAP(makerR1.all, '1', params.outdirAnnotate)
        .set { snap1_ch }

    // Sort output
    // annotating.flatten().filter( ~/.*fasta/ )
    //     .map { it -> [ it.baseName[2..-1], it ]}
    //     .mix(r1_evidence_ch).mix(snap1_ch).groupTuple()
    //     .flatten().branch {
    //         name: !(it =~ /\./)
    //         chr: it =~ /fasta/
    //         est: it =~ /est2genome/
    //         protein: it =~ /protein2genome/
    //         snap: it =~ /snap/
    //         repeats: it =~ /repeat/
    //     }.set { r2_ch } // Looks like you figured out how to filter them easily...

    // Second round
    // MAKER_R2(r2_ch.name, r2_ch.chr, r2_ch.est, r2_ch.protein, r2_ch.snap, r2_ch.repeats,
    //     params.makerR2)
    //     .set { makerR2 }
    // round2_out = GET_GFF2(makerR2, '2', params.outdirAnnotate)

}
