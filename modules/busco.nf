process BUSCO {
    publishDir "$outdir/${assembler}", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/busco'

    input:
    tuple val(assembler), path(assembly)
    val(outdir)
    //
    output:
    path("*_$assembler")
    //
    exec:
    def run = assembly.baseName.replaceAll(/_.*/, '') + "_$assembler"

    script:
    """
    busco -i $assembly \
    -l sordariomycetes_odb10 \
    -o $run\
    -m genome \
    --offline \
    --download_path $projectDir/busco_downloads
    """
}

process EXTRACT_BUSCO {
    input:
    path("*")
    //
    output:

    //
    script:
    """
    cd
    """
    //
}
