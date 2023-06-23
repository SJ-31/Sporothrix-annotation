/*
 * Module imports
 */

include { REPEATMODELER; GET_REPEAT } from '../modules/repeatmodeler'

/*
 * Main workflow
 */

workflow repeats {
    take:
    genome_ch

    main:
    rm_ch = REPEATMODELER(genome_ch)
    GET_REPEAT(rm_ch, params.repeats, params.repeat_ref)

}
