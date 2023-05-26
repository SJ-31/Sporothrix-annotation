process BUSCO {
    publishDir "$outdir"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/busco'

    input:
    tuple val(name), path(assembly)
    val(mode)
    val(outdir)
    //
    output:
    path("${assembly.baseName.replaceAll(/_.*/, '')}")
    //
    exec:
    def run = assembly.baseName.replaceAll(/_.*/, '')

    script:
    """
    busco -i $assembly \
    -l sordariomycetes_odb10 \
    -o $run\
    -m $mode \
    --offline \
    --download_path $projectDir/busco_downloads
    """
}

process EXTRACT_BUSCO {
    publishDir "$outdir", mode: 'copy'

    input:
    path("*")
    val(outdir)
    //
    output:
    path("short_summary.*.txt")
    //
    shell:
    '''
    cp */short_summary.*.txt .
    '''
    //
}
