process ONTOLOGIZER {
    publishDir "$outdir", mode: "copy"

    input:
    path(sample_stats)
    val(outdir)
    //

    output:
    path("${sample_stats.baseName}_GO_enrichment.txt")
    //

    script:
    """
    Rscript ontologizer_prep.r $sample_stats $params.go_annotations_ref
    java -jar $params.ontologizer_jar -g $params.go_ontology \
        -a ontologizer_mappings.ids \
        -p ${sample_stats.baseName}_population.txt \
        -s ${sample_stats.baseName}_interest.txt \
        -m "Bonferroni-Holm"
    mv *table* ${sample_stats.baseName}_GO_enrichment.txt
    """
    //
}
