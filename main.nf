/*
 * Raw files for assembly, if you don't split up the dataset, will run out of memory
 */
// params.raw = "$projectDir/data/raw1-4" // Completed Sat 13 May, 2023
// params.raw = "$projectDir/data/raw5-7" // Completed Sat 13 May, 2023
// params.raw = "$projectDir/data/raw8" // Completed Sun 14 May, 2023
// params.raw = "$projectDir/data/raw9-10" // Completed Sat 13 May, 2023
// params.raw = "$projectDir/data/raw11-12" // Completed Sun 15 May, 2023
// params.raw = "$projectDir/data/raw13-16" // Completed Sat 13 May, 2023

params.results = "$projectDir/results"
// raw_ch = Channel.fromFilePairs("$params.raw/S*_R{1,2}_001.fastq.gz")

/*
 * Assembly file channels
 */
megahit = "$params.results/assembly/3-megahit/*"
spades = "$params.results/assembly/3-spades/*"
testing = "$params.results/assembly/test_annotate/*"

assembly_ch = Channel.fromPath(
    // "$megahit",
    // "$spades",
    "$testing")
    .map {it -> [ it.baseName[2..-1], it ]} // Remove the prefix for sorting purposes

/*
 * Workflow imports
 */

include { assembly } from './workflows/assembly'
include { annotation } from './workflows/annotation'

/*
 *
 */
workflow {
    // assembly(raw_ch) Completed for now, Mon 15 May, 2023
    annotation(assembly_ch)
}
