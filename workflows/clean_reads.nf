include { FASTQC as FASTQC_I } from '../modules/fastqc'
include { FASTQC as FASTQC_T } from '../modules/fastqc'
include { FASTP } from '../modules/fastp'
include { BBDUK } from '../modules/bbduk'
include { COMBINE as COMBINE_QC } from '../modules/combine'

workflow clean_reads {
    take:
    raw

    main:
    FASTQC_I(raw, "$params.assembly/1-initialChecks")
    // Initial quality check
        .zip.collect().set { fastqcInitial_ch }
    COMBINE_QC(fastqcInitial_ch)
    FASTP(raw, "$params.assembly/2-postTrim") // Check
        .set { fastp }
    FASTQC_T(fastp.fastq, params.outdirTrim) // Verify results of trimming
        .zip.collect().set { fastqcTrim_ch }
    BBDUK(fastp.fastq, params.mtDNA,
    params.bbduk_args,
    "$params.assembly/2.5-bbduk_flagged_seqs",
    params.cleaned)
        .set { bbduk_ch }
}
