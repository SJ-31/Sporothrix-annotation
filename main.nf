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
 * Assembly file channels
 */
current = "$params.results/assembly/annotate-current/*"
to_annotate = "$params.results/assembly/annotate/*"

assembly_ch = Channel.fromPath(
    "$current")
    .map {it -> [ it.baseName[2..-1], it ]} // Add a prefix for sorting purposes

assess_ch = Channel.fromPath(
    "$to_annotate")
    .map {it -> [ it.baseName[2..-1], it ]} // Add a prefix for sorting purposes

/*
 * Workflow imports
 */

include { assembly; assess } from './workflows/assembly'
include { annotation } from './workflows/annotation'
include { rnaseq } from './workflows/rnaseq'

/*
 *
 */
workflow {
    // assembly(raw_ch) Completed for now, Mon 15 May, 2023
    // rnaseq(rna_ch) Completed Fri 19 May, 2023
    // annotation(assembly_ch)
    // assess(assess_ch)
}
