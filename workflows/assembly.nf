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
include { BUSCO; EXTRACT_BUSCO; MULTIQC_B } from '../modules/busco'

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
    EXTRACT_BUSCO(
        BUSCO(spadesA_ch.mix(megahitA_ch), params.buscoOut)
    ).branch{
        spades: it =~ /spades/
        megahit: it =~ /megahit/
    }.set { busco_ch }
    MULTIQC_B(busco_ch.spades.collect()
            .mix(busco_ch.megahit.collect()),
            params.mqcOutdirA)
    // MULTIQC_B(busco_ch, params.mqcOutdirA)


    emit:
    spadesA_ch
    megahitA_ch
}
