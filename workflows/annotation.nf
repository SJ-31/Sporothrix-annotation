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

workflow annotation {
    take:
    assemblies

    main:
    genemark_ch = GENEMARKS_ES(assemblies)
    assemblies.mix(genemark_ch).groupTuple(sort: true)
        .set { r1_ch }
    MAKER_R1(r1_ch, params.makerR1, true)
        .set { makerR1 }
    round1_out = GET_GFF(makerR1, '1')
    round1_out.evidence.transpose()
        .set { r1_evidence_ch }

    // Train predictors
    snap1_ch = SNAP(round1_out.all, '1', params.outdirAnnotate).hmm
    assemblies.mix(round1_out.all).groupTuple(sort: true)
        .set { to_aug1 }
    aug1_ch = AUGUSTUS(to_aug1, '1', params.outdirAnnotate) // Will use the S. schenckii reference genome for augustus training
    assemblies.mix(r1_evidence_ch).mix(snap1_ch)
        .groupTuple(sort: true)
        .set { r2_ch }

    // Second round
    // MAKER_R2(r2_ch, params.makerR2, false)

}
