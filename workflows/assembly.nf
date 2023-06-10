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
include { BUSCO; EXTRACT_BUSCO } from '../modules/busco'
include { QUAST } from '../modules/quast'

/*
 *
 */

workflow assembly {
    take:
    cleaned

    main:
    if  ( params.withSpades ) {
        EXTRACT_SPADES(SPADES(cleaned, params.spadesargs), params.spadesOut)
            .set{ spades_ch }
    } else { spades_ch = Channel.empty()
    }
    if ( params.withMegahit ) {
        EXTRACT_MH(MEGAHIT(cleaned), params.megaOut)
            .set { megahit_ch }
    } else { megahit_ch = Channel.empty()
    }
    assemblies = megahit_ch.mix(spades_ch)

    // Quality assessment
    if ( params.assess_assemblies )
        BUSCO(assemblies, 'genome', params.buscoOutA)
        .branch{
            spades: it =~ /spades/
            megahit: it =~ /megahit/
        }.set { busco_ch }
        // assemblies.flatten().branch {
        //     fasta: it =~/fasta/
        // }.set { genomes }
        // QUAST(genomes.fasta.collect(), params.quastOut, params.quastRef,
        // params.quastRefF, params.quast_args)
}
