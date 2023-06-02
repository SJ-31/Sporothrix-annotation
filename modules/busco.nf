process BUSCO {
    publishDir "$outdir/${assembler}", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/busco'

    input:
    tuple val(assembler), path(assembly)
    val(mode)
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
    -m $mode \
    --offline \
    --download_path $projectDir/busco_downloads
    """
}

process EXTRACT_BUSCO {
    input:
    path("*")
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

process MAKER_BUSCO {
    publishDir "$outdir/${sample}", mode: 'copy'
    conda '/home/sc31/Bio_SDD/miniconda3/envs/busco'

    input:
    tuple val(sample), path(assembly)
    val(mode)
    val(outdir)
    //
    output:
    path("*_BUSCO")
    //
    script:
    """
    busco -i $assembly \
    -l sordariomycetes_odb10 \
    -o ${assembly}_BUSCO \
    -m $mode \
    --offline \
    --download_path $projectDir/busco_downloads
    """
}

process BUSCO_F {
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

process EXTRACT_BUSCO_F {
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
