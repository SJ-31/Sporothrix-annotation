/*
 * Workflow imports
 */

include { assembly; assess } from './workflows/assembly'
include { train_genemarks; annotation; get_buscos } from './workflows/annotation'
include { repeats } from './workflows/repeatlibrary'
include { rnaseq } from './workflows/rnaseq'
include { scaffold; assess_scaffolds } from './workflows/finishing'

/*
 * Raw files for assembly, if you don't split up the dataset, will run out of memory
 */
params.raw = "$projectDir/data/$params.to_assemble/"
params.rna = "$projectDir/data/raw_rna"
params.results = "$projectDir/results"

raw_ch = Channel.fromFilePairs("$params.raw/S*_R{1,2}_001.fastq.gz")
rna_ch = Channel.fromFilePairs("$params.rna/*_{1,2}.fastq")

/*
 * Repeat library channels
 */

Channel.fromPath("$projectDir/data/reference/genomes/for_repeats/*")
    .map { it -> [ it.baseName, it ]}
    .set { genome_ch }

/*
 * Finishing channels
 */

Channel.fromPath("$projectDir/results/assembly/$params.to_scaffold/*")
    .map { it -> [ it.baseName.replaceAll(/_.*/, '').replaceAll(/.*-/, '')
                , it ]}
    .set { contigs_ch }
contigs_ch.combine(raw_ch, by: 0)
    .set { contigs_reads_ch }
Channel.fromPath("$projectDir/data/reference/genomes/scaffold_ref/*")
    .collect().set { ref_ch }
Channel.fromPath(params.genomes)
    .map { it -> it.baseName + '.fasta,2' }
    .collectFile(name: 'ref.txt', newLine: true )
    .set { ref_config }

/*
 * Annotation channels
 */
annotation = "$params.results/assembly/7-aligned/$params.current/*"
Channel.fromPath("$params.results/assembly/5-scaffolds/ragout/*.fasta")
    .map { it -> [ it.baseName, it ] }
    .set { all_scaffs }
Channel.fromPath(
    "$annotation")
    .map {it -> [ it.baseName, it ]}
    .set { chromosome_ch  }
Channel.fromPath(
    "$params.results/assembly/5-scaffolds/ragout/${params.current}_scaffolds.fasta")
    .map { it -> [ it.baseName, it ] }
    .set { scaffold_ch }

Channel.fromPath(
    "$params.results/annotation/S*/S*BUSCO")
    .map { it -> [ it.baseName.replaceAll(/-.*/), it ]}
    .set { busco_results_ch }
/*
 *
 */
workflow {
    if ( params.assemble_genome )
        assembly(raw_ch)
    if ( params.assess_assemblies )
        assess(assess_ch)
    if ( params.assemble_transcriptome )
        rnaseq(rna_ch)
    if ( params.get_replib )
        repeats(genome_ch)
    if ( params.scaffold_contigs )
        scaffold(contigs_ch, contigs_reads_ch, ref_ch )
    if ( params.assess_scaffolds )
        assess_scaffolds()
    // train_genemarks(all_scaffs)
    if ( params.annotate_scaffold )
        annotation(chromosome_ch, scaffold_ch)
    if ( params.find_buscos )
        get_buscos(busco_results_ch)
}
