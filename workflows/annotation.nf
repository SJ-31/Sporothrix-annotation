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
    MAKER_F(chromosomes, params.makerR1)
        .set { makerR1 }
    round1_out = GET_GFF(makerR1, '1', params.outdirAnnotate)
    round1_out.evidence.transpose()
        .set { r1_evidence_ch }

    // Train predictors
    snap1_ch = SNAP(round1_out.all, '1', params.outdirAnnotate).hmm
    chromosomes.flatten().filter( ~/.*fasta/ )
        .map { it -> [ it.baseName[2..-1], it ]}
        .mix(r1_evidence_ch).mix(snap1_ch)
        .groupTuple(sort: true)
        .set { r2_ch }

    // Second round
    MAKER_R2(r2_ch, params.makerR2)
        .set { makerR2 }
    round2_out = GET_GFF2(makerR2, '2', params.outdirAnnotate)

}
