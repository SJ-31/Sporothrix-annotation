process ADD_RG {
    input:
    tuple val(pair_id), path(aligned_sam)
    val(params)
    //
    output:
    tuple val(pair_id), path("*RG*")
    //
    script:
    """
    java -jar /home/sc31/Bio_SDD/tools/picard.jar AddOrReplaceReadGroups \
    I=$aligned_sam \
    O=${pair_id}_aligned_withRG.sam \
    $params
    """
    //
}
