/*
 * Module imports
 */
include { BUSCO } from "../modules/busco"
include { NTJOIN } from "../modules/ntjoin"
include { SSPACE } from "../modules/sspace"
include { EXTRACTCHR } from "../modules/extractchr"
include { RAGOUT; EXTRACT_RAGOUT } from "../modules/ragout"
include { QUAST } from "../modules/quast"
include { CAP } from "../modules/cap"
include { FILTER_MT } from "../modules/filtermt"

/*
 * Workflow
 */

workflow scaffold {
    take:
    contigs
    references

    main:
    EXTRACT_RAGOUT(
    RAGOUT(contigs, references, "$params.scaffolds"),
    params.scaffolds)
        .set { extracted_ch }
    FILTER_MT(extracted_ch, params.mtDNA, "$params.assembly/6.5-mtDNA")
    .filtered // A precaution in case BBduk filtering missed some reads, but is usually empty
        .set { filtered_ch }
    EXTRACTCHR(filtered_ch, params.quastRef, "$params.assembly/7-aligned")
}

workflow assess_scaffolds {
    extracted_ch = Channel.fromPath("$params.scaffolds/ragout/*.fasta")
    BUSCO(extracted_ch.map { it -> [ it.baseName.replaceAll(/.*-/, ''), it ]},
    'genome', "$params.assembly/6-assess_scaff")
}
