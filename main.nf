/*
 * Workflow imports
 */

include { assembly; assess } from './workflows/assembly'
include { annotation } from './workflows/annotation'
include { repeats } from './workflows/repeatlibrary'
include { rnaseq } from './workflows/rnaseq'
include { scaffold; assess_scaffolds } from './workflows/finishing'

/*
 * Raw files for assembly, if you don't split up the dataset, will run out of memory
 */
params.raw = "$projectDir/data/$params.to_assemble"
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
 * Assembly file channels
 */
current = "$params.results/assembly/annotate-current/*"
to_annotate = "$params.results/assembly/annotate/*"

Channel.fromPath(
    "$current")
    .map {it -> [ it.baseName[2..-1], it ]}
    .set { assembly_ch } // Add a prefix for sorting purposes

assess_ch = Channel.fromPath(
    "$to_annotate")
    .map {it -> [ it.baseName[2..-1], it ]} // Add a prefix for sorting purposes

/*
 * Finishing channels
 */

Channel.fromPath("$projectDir/results/assembly/annotate/*")
    .map { it -> [ (it =~ /.*-(.*)_.*/)[0][1], it ]}
    .set { contigs_ch }
Channel.fromPath("$projectDir/results/assembly/annotate-current/*")
    .map { it -> [ (it =~ /.*-(.*)_.*/)[0][1], it ]}
    .set { test_ch }
contigs_ch.combine(raw_ch, by: 0)
    .set { contigs_reads_ch }
Channel.fromPath("$projectDir/data/reference/genomes/scaffold_ref/*")
    .collect().set { ref_ch }
Channel.fromPath(params.genomes)
    .map { it -> it.baseName + '.fasta,2' }
    .collectFile(name: 'ref.txt', newLine: true )
    .set { ref_config }
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
        // scaffold(test_ch, contigs_reads_ch, ref_ch ) // For testing purposes
    if ( params.assess_scaffolds )
        assess_scaffolds()
    if ( params.annotate_scaffold )
        annotation(assembly_ch)
}
