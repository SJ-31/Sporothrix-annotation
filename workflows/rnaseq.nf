/*
 * Module imports
 */
include { FASTQC as FASTQC_I } from '../modules/fastqc'
include { FASTQC as FASTQC_T } from '../modules/fastqc'
include { MULTIQC as MULTIQC_I } from '../modules/multiqc'
include { MULTIQC as MULTIQC_T } from '../modules/multiqc'
include { COMBINE as COMBINE_QC } from '../modules/combine'
include { FASTP } from '../modules/fastp'
include { RNA_SPADES; EXTRACT_RNASPADES } from '../modules/spades'
include { BUSCO } from '../modules/busco'
include { QUAST } from '../modules/quast'

/*
 * Main workflow
 */

workflow rnaseq {
    take:
    raw

    main:
    FASTQC_I(raw, "$params.rnaseq/1-initialChecks") // Initial quality check
        .zip.collect().set { fastqcInitial_ch }
    COMBINE_QC(fastqcInitial_ch)
    FASTP(raw, "$params.rnaseq/2-postTrim") // Check
        .set { fastp }
    FASTQC_T(fastp.fastq, "$params.rnaseq/2-postTrim") // Verify results of trimming
        .zip.collect().set { fastqcTrim_ch }
    spades_ch = RNA_SPADES(fastp.fastq)
    EXTRACT_RNASPADES(spades_ch, "$params.rnaseq/3-spades")
        .set{ spadesA_ch }
    // Quality assessment
    BUSCO(spadesA_ch, 'transcriptome', "$params.rnaseq/4-assess")
}
