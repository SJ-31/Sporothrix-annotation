/*
 * Module imports
 */
include { FASTQC as FASTQC_I } from '../modules/fastqc'
include { FASTQC as FASTQC_T } from '../modules/fastqc'
include { MULTIQC as MULTIQC_I } from '../modules/multiqc'
include { MULTIQC as MULTIQC_T } from '../modules/multiqc'
include { COMBINE as COMBINE_QC } from '../modules/combine'
include { FASTP } from '../modules/fastp'
include { COMBINE as COMBINE_P } from '../modules/combine'
include { SEQKIT_STATS } from '../modules/seqkit'
include { SPADES } from '../modules/spades'

/*
 *
 */

workflow assembly {
    take:
    raw

    main:
    FASTQC_I(raw, params.outdirInitial) // Initial quality check
        .zip.collect().set { fastqcInitial_ch }
    MULTIQC_I(COMBINE_QC(fastqcInitial_ch), params.mqcOutdirI)
    FASTP(raw, params.outdirTrim) // Check
        .set { fastp }
    FASTQC_T(fastp.fastq, params.outdirTrim) // Verify results of trimming
        .zip.collect().set { fastqcTrim_ch }
    MULTIQC_T(COMBINE_P(fastp.html.mix(fastp.json)
                        .mix(fastqcTrim_ch)
                        .collect()), params.mqcOutdirT)
    SPADES(fastp.fastq, params.spadesOut).view()

}
