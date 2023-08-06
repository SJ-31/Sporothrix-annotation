/*
 * Module imports
 */

include { SPADES; EXTRACT_SPADES } from '../modules/spades'
include { MEGAHIT; EXTRACT_MH } from '../modules/megahit'
include { BUSCO } from '../modules/busco'
include { QUAST } from '../modules/quast'

workflow assembly {
    take:
    cleaned

    main:
    if  ( params.withSpades ) {
        EXTRACT_SPADES(SPADES(cleaned, params.spadesargs),
        "$params.assembly/3-spades")
            .set{ spades_ch }
    } else { spades_ch = Channel.empty()
    }
    if ( params.withMegahit ) {
        EXTRACT_MH(MEGAHIT(cleaned),
        "$params.assembly/3-megahit")
            .set { megahit_ch }
    } else { megahit_ch = Channel.empty()
    }
    assemblies = megahit_ch.mix(spades_ch)

    // Quality assessment
    if ( params.assess_assemblies )
        BUSCO(assemblies, 'genome',
        "$params.assemblies/4-assess_assemblies")
        .branch{
            spades: it =~ /spades/
            megahit: it =~ /megahit/
        }.set { busco_ch }
        assemblies.flatten().branch {
            fasta: it =~/fasta/
        }.set { genomes }
        QUAST(genomes.fasta.collect(), "$params.assemblies/4-assess_assemblies",
        params.reference,
        params.ref_gff, params.quast_args)
}
