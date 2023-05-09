/*
 * Params
 */
params.raw = "$projectDir/data/raw"
params.results = "$projectDir/results"


raw_ch = Channel.fromFilePairs("$params.raw/S*_R{1,2}_001.fastq.gz")

/*
 * Workflow imports
 */

include { assembly } from './workflows/assembly'

/*
 *
 */
workflow {
    assembly(raw_ch)

}
