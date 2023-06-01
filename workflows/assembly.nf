/*
 * Module imports
 */
include { FASTQC as FASTQC_I } from '../modules/fastqc'
include { FASTQC as FASTQC_T } from '../modules/fastqc'
include { MULTIQC as MULTIQC_I } from '../modules/multiqc'
include { MULTIQC as MULTIQC_T } from '../modules/multiqc'
include { COMBINE as COMBINE_QC } from '../modules/combine'
include { FASTP } from '../modules/fastp'
include { BBDUK } from '../modules/bbduk'
include { COMBINE as COMBINE_P } from '../modules/combine'
include { SEQKIT_STATS } from '../modules/seqkit'
include { SPADES; EXTRACT_SPADES } from '../modules/spades'
include { MEGAHIT; EXTRACT_MH } from '../modules/megahit'
include { BUSCO; EXTRACT_BUSCO; MULTIQC_B } from '../modules/busco'
include { QUAST } from '../modules/quast'

/*
 *
 */

workflow assembly {
    take:
    raw

    main:
    // FASTQC_I(raw, params.outdirInitial) // Initial quality check
    //     .zip.collect().set { fastqcInitial_ch }
    // COMBINE_QC(fastqcInitial_ch)
    FASTP(raw, params.outdirTrim) // Check
        .set { fastp }
    BBDUK(fastp.fastq, params.mtDNA,
            params.bbduk_args, params.bbdukOut).set { bbduk_ch }
    // FASTQC_T(fastp.fastq, params.outdirTrim) // Verify results of trimming
    //     .zip.collect().set { fastqcTrim_ch }
    if  ( params.withSpades ) {
        spades_ch = SPADES(fastp.fastq, params.spadesargs)
        EXTRACT_SPADES(spades_ch, params.spadesOut)
            .set{ spadesA_ch }
    }
    if ( params.withMegahit ) {
        megahit_ch = MEGAHIT(bbduk_ch.reads)
        EXTRACT_MH(megahit_ch, params.megaOut)
    }
    // Quality assessment
}

workflow assess {
    take:
    assemblies

    main:
    EXTRACT_BUSCO(
        BUSCO(assemblies, 'genome', params.buscoOut)
    ).branch{
        spades: it =~ /spades/
        megahit: it =~ /megahit/
    }.set { busco_ch }

    assemblies.flatten().branch {
        fasta: it =~/fasta/
    }.set { genomes }

    QUAST(genomes.fasta.collect(), params.quastOut, params.quastRef, params.quastRefF,
        params.quast_args)
}
