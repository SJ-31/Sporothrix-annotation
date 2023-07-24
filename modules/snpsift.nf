process SNPSIFT {
    publishDir "$outdir/variants/$pair_id", mode: 'copy'

    input:
    tuple val(pair_id), path(annotations)
    val(outdir)
    //
    output:
    tuple (val(pair_id), path("${pair_id}_filtered_GO.vcf"))
    //
    script:
    """
    snpsift_wrapper.sh $pair_id $annotations $params.snpsift_jar
    Rscript $bin/go_anno.r ${pair_id}_annotated_vars.vcf $params.go_annotations_ref
    """
    //
}
