/*
 * Workflow imports
 */
include { assembly } from './workflows/assembly'
include { train_genemarks; annotation; get_buscos } from './workflows/annotation'
include { repeats } from './workflows/repeatlibrary'
include { variant_calling; extract_buscos } from './workflows/variant_calling'
include { rnaseq } from './workflows/rnaseq'
include { scaffold; assess_scaffolds } from './workflows/finishing'
include { clean_reads } from './workflows/clean_reads.nf'
/*
 * Raw files for assembly, if you don't split up the dataset, will run out of memory
 */
raw_ch = Channel.fromFilePairs("$params.raw/S*_R{1,2}_001.fastq.gz")
rna_ch = Channel.fromFilePairs("$params.rna/*_{1,2}.fastq")
clean_ch = Channel.fromFilePairs("$params.clean/B-S*_R{1,2}_001.fastq.gz")

/*
 * Repeat library channels
 */
Channel.fromPath("$projectDir/data/reference/genomes/for_repeats/*")
    .map { it -> [ it.baseName, it ]}
    .set { genome_ch }

/*
 * Finishing channels
 */
Channel.fromPath("$projectDir/results/assembly/$params.to_scaffold")
    .map { it -> [ it.baseName.replaceAll(/_.*/, '').replaceAll(/.*-/, ''), it ]}
    .set { contigs_ch }
Channel.fromPath("$projectDir/data/reference/genomes/scaffold_ref/*")
    .collect().set { ref_ch }

/*
 * Annotation channels
 */
annotation = "$params.assembly/7-aligned/$params.current/*"
Channel.fromPath("$params.scaffolds/ragout/*.fasta")
    .map { it -> [ it.baseName, it ] }
    .set { train_ch }
Channel.fromPath(
    "$annotation")
    .map {it -> [ it.baseName, it ]}
    .set { chromosome_ch  }
Channel.fromPath(
    "$params.scaffolds/${params.current}_scaffolds.fasta")
    .map { it -> [ it.baseName, it ] }
    .set { scaffold_ch }
Channel.fromPath(
    "$params.annotation/S*/S*BUSCO")
    .map { it -> [ it.baseName.replaceAll(/-.*/), it ]}
    .set { busco_results_ch }

/*
 * Variant calling channels
 */
Channel.fromPath(params.vc_ref)
    .set { vc_ref_ch }
Channel.fromFilePairs("$projectDir/data/cleaned_dna/${params.called_reads}/B-S*_R{1,2}_001.fastq.gz")
    .set { call_reads_ch }

workflow {
    /*
     * Assembly
     */
    if ( params.clean_reads )
        clean_reads(raw_ch)
    if ( params.assemble_genome )
        assembly(clean_ch)
    if ( params.assemble_transcriptome )
        rnaseq(rna_ch)
    if ( params.get_replib )
        repeats(genome_ch)
    if ( params.scaffold_contigs )
        scaffold(contigs_ch, ref_ch )
    if ( params.assess_scaffolds )
        assess_scaffolds()

    /*
     * Annotation
     */
    if ( params.train_genemarks )
        train_genemarks(train_ch)
    if ( params.annotate_scaffold )
        annotation(chromosome_ch, scaffold_ch)
    if ( params.collect_buscos )
    // Combine busco gene sequences from busco runs on different sample
        get_buscos(busco_results_ch)

    /*
     * Variant calling
     */
    if ( params.call )
        variant_calling(vc_ref_ch, call_reads_ch)
    if ( params.extract_busco_genes )
    // Extract busco sequences from aligned BAM file and generate a multiple sequence
    //  alignment
        extract_buscos(params.vc_ref, params.vc_align, "$params.cleaned/*")
}
