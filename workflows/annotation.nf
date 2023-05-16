/*
 * Module imports
 */
include { GET_GFF } from "../modules/maker"
include { MAKER as MAKER_R1 } from "../modules/maker"
include { MAKER as MAKER_R2 } from "../modules/maker"
include { MAKER as MAKER_R3 } from "../modules/maker"
include { SNAP } from "../modules/snap"
include { GENEMARKS_ES } from "../modules/genemarksES"

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
    makerR1.view()
    // round1_out = GET_GFF(makerR1)
    // SNAP(round1_out.all)
    // round1.view()

}
