process LIFTOFF {
    publishDir "$outdir", mode: 'copy', pattern: "*gff"
    conda '/home/sc31/Bio_SDD/miniconda3/envs/liftoff'

    input:
    tuple val(name), path(target)
    path(reference)
    path(reference_gff)
    val(outdir)
    //
    output:
    tuple val(name), path("${name}_lifted.gff"), path(reference)
    //
    script:
    """
    liftoff $target $reference -g $reference_gff -o ${name}_lifted.gff
    """
    //
}