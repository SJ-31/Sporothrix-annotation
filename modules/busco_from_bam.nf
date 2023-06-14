process BUSCO_FROM_BAM {
    input:
    tuple val(name), path(sample)
    path(reference)
    val(busco_ref)
    //
    output:
    path("*at*fasta")
    //
    shell:
    '''
    minimap2 -a !{reference} !{sample} > aligned.sam
    while read line
        do
        fields=( $(echo $line | awk '! /#/ && $2 ~ /Complete/ {print $1,$3,$4,$5,$6,$10}') )
        case ${fields[4]} in
        '-')
            region=${fields[1]}:${fields[3]}-${fields[2]};;
        '+')
            region=${fields[1]}:${fields[2]}-${fields[3]};;
        *)
            continue;;
        esac
        echo ">!{name}_${fields[0]}:${fields[5]}" > \
        !{name}_${fields[0]}.fasta
        samtools view aligned.bam $region | \
        awk '{print $10}'>> !{name}_${fields[0]}.fasta
        done < !{busco_ref}
    '''
    // fields[0] = Gene busco ID
    // fields[1] = Gene sequence location
    // fields[2], fields[3] = Gene start, stop
    // fields[4] = Gene strandedness
    // fields[5] = Gene description
}