process BUSCO {
    publishDir "$outdir/${assembler.baseName}", mode: 'copy',
    // conda '/home/sc31/Bio_SDD/miniconda3/envs/busco'

    input:
    tuple val(assembler), val(run), path(assembly)
    val(outdir)
    //
    output:
    path(run)
    //
    script:
    """
    busco -i $assembly
    -l sordariomycetes_odb10 \
    -o $run\
    -m genome \
    --offline \
    --download_path $projectDir/busco_downloads
    """
    //
}
