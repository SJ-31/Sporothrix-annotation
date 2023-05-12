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
include { SPADES; EXTRACT_SPADES } from '../modules/spades'
include { MEGAHIT; EXTRACT_MH } from '../modules/megahit'
include { BUSCO } from '../modules/busco'

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
    spades_ch = SPADES(fastp.fastq)
    EXTRACT_SPADES(spades_ch, params.spadesOut)
        .set{ spadesA_ch }
    megahit_ch = MEGAHIT(fastp.fastq)
    EXTRACT_MH(megahit_ch, params.megaOut)
        .set { megahitA_ch }
    // Quality assessment
    // assembly_MH_ch.mix(assembly_S_ch).view()
    BUSCO(spadesA_ch.mix(megahitA_ch), params.buscoOut).view()

    emit:
    spadesA_ch
    megahitA_ch
}
