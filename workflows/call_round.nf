include { HAPLOTYPECALLER }  from '../modules/haplotypecaller'
include { FILTERSNPS }  from '../modules/filtersnps'
include { FILTERINDELS }  from '../modules/filterindels'
include { SELECTVARIANTS }  from '../modules/selectvariants'
include { FILTERVARIANTS } from '../modules/filtervariants'

workflow call_variants {
    take:
    sorted_bams
    refs
    round

    main:
    HAPLOTYPECALLER(sorted_bams, refs, round)
        .set { hc1_output_ch }
    SELECTVARIANTS(hc1_output_ch, refs, round)
        .set { sv_ch }
    FILTERSNPS(sv_ch.raw_snps_ch, refs, round, params.vc)
        .set { snps_ch }
    FILTERINDELS(sv_ch.raw_indels_ch, refs, round, params.vc)
        .set { indels_ch }
    if ( round != 'final') {
        FILTERVARIANTS(snps_ch, indels_ch, refs, round)
            .set { filtered_variants_ch }
    } else {
        filtered_variants_ch = Channel.empty()
    }
    emit:
    snps_ch
    indels_ch
    filtered_variants_ch
}
