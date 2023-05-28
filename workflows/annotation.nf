/*
 * Module imports
 */
include { GET_GFF } from "../modules/maker"
include { MAKER as MAKER_R1 } from "../modules/maker"
include { MAKER as MAKER_R2 } from "../modules/maker"
include { MAKER as MAKER_R3 } from "../modules/maker"
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
    MAKER_R1(chromosomes, params.makerR1, true)
        .set { makerR1 }
    round1_out = GET_GFF(makerR1, params.outdirAnnotate, '1')
    round1_out.evidence.transpose()
        .set { r1_evidence_ch }

    // Train predictors
    snap1_ch = SNAP(round1_out.all, '1', params.outdirAnnotate).hmm
    // chromosomes.mix(r1_evidence_ch).mix(snap1_ch)
    //     .groupTuple(sort: true)
    //     .view()
    //     .set { r2_ch }

    // Second round
    // MAKER_R2(r2_ch, params.makerR2, false)

}
