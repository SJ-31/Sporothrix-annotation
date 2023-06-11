process BUSCO {
    publishDir "$outdir/${sample}", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/busco'

    input:
    tuple val(sample), path(fasta)
    val(mode)
    val(outdir)
    //
    output:
    path("*_BUSCO")
    //
    exec:
    def name = fasta.baseName

    script:
    """
    busco -i $fasta \
    -l sordariomycetes_odb10 \
    -o ${name}_BUSCO \
    -m $mode \
    --offline \
    --download_path $projectDir/busco_downloads
    """
}
