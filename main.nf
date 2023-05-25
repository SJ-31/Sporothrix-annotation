/*
 * Raw files for assembly, if you don't split up the dataset, will run out of memory
 */
// params.raw = "$projectDir/data/raw1-4" // Completed Sat 13 May, 2023
// params.raw = "$projectDir/data/raw5-7" // Completed Sat 13 May, 2023
// params.raw = "$projectDir/data/raw8" // Completed Sun 14 May, 2023
// params.raw = "$projectDir/data/raw9-10" // Completed Sat 13 May, 2023
// params.raw = "$projectDir/data/raw11-12" // Completed Sun 15 May, 2023
// params.raw = "$projectDir/data/raw13-16" // Completed Sat 13 May, 2023
params.rna = "$projectDir/data/raw_rna"

params.results = "$projectDir/results"
// raw_ch = Channel.fromFilePairs("$params.raw/S*_R{1,2}_001.fastq.gz")

// Fri 19 May, 2023 Assembled SRR9602168
rna_ch = Channel.fromFilePairs("$params.rna/*_{1,2}.fastq")

/*
 * Workflow imports
 */

include { assembly; assess } from './workflows/assembly'
include { annotation } from './workflows/annotation'
include { repeats } from './workflows/repeatlibrary'
include { rnaseq } from './workflows/rnaseq'
include { scaffold } from './workflows/finishing'
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

Channel.fromPath("$projectDir/results/annotate/*")
    .map { it -> [ (it =~ /.*-(.*)_.*/)[0][1], it ]}
    .set { contigs_ch }
contigs_ch.combine(raw_ch, by: 0)
    .set { contigs_reads_ch }

/*
 *
 */
workflow {
    // assembly(raw_ch) // Completed for now, Mon 15 May, 2023
    // rnaseq(rna_ch) // Completed Fri 19 May, 2023
    // repeats() Completed Sun 21 May, 2023
    scaffold(contigs_ch, contigs_reads_ch)
    // annotation(assembly_ch)
    // assess(assess_ch)

}
