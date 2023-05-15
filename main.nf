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

/*
 * Fasta file channels
 */


raw_ch = Channel.fromFilePairs("$params.raw/S*_R{1,2}_001.fastq.gz")
// results_ch = Channel.fromPath("$params.")
// results_ch.view()
raw_ch.view()

/*
 * Workflow imports
 */

include { assembly } from './workflows/assembly'
include { }

/*
 *
 */
workflow {
    // assembly(raw_ch)
}
