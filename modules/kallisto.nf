process KALLISTO {
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(sample), path(combined), path(reads)
    val(outdir)
    //
    output:
    tuple val(sample), path("${sample}_kallisto")
    //
    script:
    """
    kallisto index -i index $combined
    kallisto quant -i index -o ${sample}_kallisto ${reads[0]} ${reads[1]}
    """
    //
}
