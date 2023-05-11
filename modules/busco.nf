process BUSCO {
    publishDir "$outdir", mode: 'copy',
    conda '/home/sc31/Bio_SDD/miniconda3/envs/busco'

    input:
    //
    output:
    //
    script:
    """
    busco -i
    -l sordariomycetes_odb10 \
    -o \
    -m \
    --offline \
    --download_path $projectDir/busco_downloads
    """
    //
}
