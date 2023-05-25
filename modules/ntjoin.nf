process NTJOIN {
    conda "/home/sc31/Bio_SDD/miniconda3/envs/ntjoin"
    publishDir "$outdir", mode: 'copy'
    input:
    tuple val(name), path(contigs)
    path(references)
    each path(config)
    val(outdir)
    //
    output:
    tuple val(name), path("*all.scaffolds.fa")
    //
    script:
    """
    ntJoin assemble target=$contigs \
    target_weight=1 \
    g=10 \
    k=25 \
    reference_config=$config
    """
    //
}
