/*  GATK4 Variant Calling Pipeline
 *  Author: Mohammed Khalfan < mkhalfan@nyu.edu >
 *  NYU Center for Genetics and System Biology 2020
 */

// Setting some defaults here,
// can be overridden in config or via command line

/*
 * Module imports
 */
include { COMBINE } from '../modules/combine'
include { CONS } from '../modules/emboss_cons'
include { MAKE_REF } from '../modules/makeref'
include { MERGE_VARIANTS } from '../modules/merge_variants'
include { ALIGN }  from '../modules/align.nf'
include { MARKDUPLICATESSPARK }  from '../modules/markduplicates'
include { ADD_RG } from '../modules/add_readgroups'
include { GETMETRICS }  from '../modules/getmetrics'
include { ANALYZECOVARIATES }  from '../modules/analyzecovariates'
include { QC }  from '../modules/qc'
include { SNPEFF }  from '../modules/snpeff'
include { BQSR } from '../modules/bqsr'
include { VCF_GET_REGION } from '../modules/region_from_vcf'
include { BUSCO_FROM_BAM } from '../modules/busco_from_bam'
include { LIFTOFF } from '../modules/liftoff'
include { MUSCLE } from '../modules/muscle'
include { NJ } from '../modules/nj'
include { MAFFT } from '../modules/mafft'
include { GFF_BUSCO_MAP } from '../modules/busco_gff_mapping'
include { KALLISTO } from '../modules/kallisto'
include { BUSCO_FROM_GFF } from '../modules/busco_from_gff'
include { call_variants as call_variants_R1 } from './call_round'
include { call_variants as call_variants_R2 } from './call_round'

workflow variant_calling {
    take:
    refs
    cleaned_pairs_ch

    main:
    MAKE_REF(refs)
        .set { ref_ch }
    ADD_RG(ALIGN(cleaned_pairs_ch, ref_ch), params.readgroups)
        .set { aligned_reads_ch }
    if ( params.caller = "GATK") {
        MARKDUPLICATESSPARK(aligned_reads_ch, params.vc)
            .set { mds } //
        GETMETRICS(mds.sorted_bam, ref_ch, params.vc)
            .set { metrics_qc_ch }

        // First round
        call_variants_R1(mds.sorted_bam, ref_ch, 'first')
            .set { round1 }
        BQSR(mds.sorted_bam.join(round1.filtered_variants_ch), ref_ch)
            .set { bqsr_ch }

        // Final round, after recalibration
        call_variants_R2(bqsr_ch.recalibrated_bams, ref_ch, 'final')
            .set { round2 }
        MERGE_VARIANTS(round2.snps_ch, round2.indels_ch, params.vc)
            .set { merged }
        VCF_GET_REGION(merged, "$params.vc/gene_vars", params.gene_locs)
        QC(mds.dedup_qc_ch.join (metrics_qc_ch).join (round1.snps_ch)
        .join (round2.snps_ch),
           params.vc)
    } else if ( params.caller = "freeBayes" ) {
        FREEBAYES(aligned_reads_ch, ref_path, params.vc)
            .set { merged }
    }
    SNPEFF(merged, "$params.vc/snpeff")
    // ANALYZECOVARIATES(bqsr_ch.analyze_covariates_in_ch)
    //     .set { analyzed_covariates_ch }
    /* Process qc creates a report for each sample.
    * Below we compile these into a single report.
    */
}

workflow extract_buscos {
    take:
    ref_path
    assembly_path
    clean_path

    main:
    Channel.fromPath(assembly_path)
    .map { it -> [ it.baseName.replaceAll(/_.*/, ''), it ]}
        .set { assembly_ch }
    Channel.fromFilePairs("$clean_path/B-S*_R{1,2}_001.fastq.gz")
    .map { it -> [ it[0].replaceAll(/.*-/, ''), it[1] ]}
        .set { reads_ch }
    LIFTOFF(assembly_ch, ref_path, params.vc_refgff, "$params.vc/1-lifted_over")
        .set { lift_ch }
    GFF_BUSCO_MAP(lift_ch, params.busco_gff,
    "$params.vc/2-liftover_busco_mapping")
        .set { busco_locs_ch }
    BUSCO_FROM_GFF(busco_locs_ch)
    .tap { per_sample_ch }
    .flatten().filter( ~/.*fasta/ )
    .map { it ->  [ it.baseName.replaceAll(/.*_/, ''), it ] }
    .groupTuple()
        .set { single_copy_buscos }
    COMBINE(single_copy_buscos.mix(per_sample_ch),
    "$params.vc/3-combined_genes").exists.branch {
        per_sample: it =~ /S.*/
        all_samples: true
    }.set { combined_ch }
    combined_ch.per_sample
    .map { it -> [ it.baseName.replaceAll(/_.*/, ''), it]}
        .set { per_sample }
    MAFFT(combined_ch.all_samples, "$params.vc/4-busco_genes_MSA")
        .set { msa_ch }
    NJ(msa_ch, params.reftree, "$params.vc/5-neighbor_joining")
    CONS(msa_ch, "$params.vc/5-consensus")
    KALLISTO(per_sample.join(reads_ch), "$params.vc/5-quantification")
}
