process FASTP {
    publishDir "$outdir/$run", mode: 'copy', pattern: "*html"
    publishDir "$outdir/$run", mode: 'copy', pattern: "*json"

    input:
    tuple val(run), path(reads)
    val(outdir)
    //
    output:
    tuple val(run), path("*.fastq.gz"), emit: fastq
    path("*.html"), emit: html
    path("*.json"), emit: json
    //
    script:
    def r1 = reads[0].baseName.replaceAll(/\..*/, '')
    def r2 = reads[1].baseName.replaceAll(/\..*/, '')
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
    --detect_adapter_for_pe \
    -c \
    -Q \
    --trim_poly_g \
    -o T-${r1}.fastq.gz -O T-${r2}.fastq.gz  \
    -z 4  \
    -h ${run}_trim.html \
    -j ${run}_fastp.json \
    -R ${run}_report
    """
    //
}
