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
    path("short_summary.*.txt")
    //
    shell:
    '''
    cp */short_summary.*.txt .
    '''
    //
}

process MULTIQC_B {
    publishDir "$outdir/$round", pattern: "*html", mode: 'copy'

    input:
    path(files)
    val(outdir)
    //
    output:
    path("*.html")
    //
    script:
    def assembler = (files[0] =~ /\d_(.*).txt/)[0][1]
    """
    mkdir $assembler
    ls --ignore=$assembler | xargs -I{} mv {} $assembler/
    multiqc $assembler
    mv *.html ${assembler}.html
    """
    //
}



