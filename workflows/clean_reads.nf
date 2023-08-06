include { FASTQC as FASTQC_I } from '../modules/fastqc'
include { FASTQC as FASTQC_T } from '../modules/fastqc'
include { FASTP } from '../modules/fastp'
include { BBDUK } from '../modules/bbduk'

workflow clean_reads {
    take:
    raw

    main:
        // Initial quality check
    FASTQC_I(raw, "$params.assembly/1-initialChecks")
    .zip.collect()
        .set { fastqcInitial_ch }
    FASTP(raw, "$params.assembly/2-postTrim") // Trim
        .set { fastp }
    FASTQC_T(fastp.fastq, "$params/2-postTrim") // Verify results of trimming
    .zip.collect()
        .set { fastqcTrim_ch }
    BBDUK(fastp.fastq, params.mtDNA, params.bbduk_args,
    "$params.assembly/2.5-bbduk_flagged_seqs", params.cleaned)
        .set { bbduk_ch }
}
