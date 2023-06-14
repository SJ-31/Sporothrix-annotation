#!/bin/bash
busco_ref=$1
gff=$2
fasta=$3
info="#BUSCO ID\tType\tSequence\tStart\tStop\tNCBI descriptor\tBUSCO descriptor"
echo -e $info > gff_busco.txt

function read_region {
    gffread -J -r "${region}" "$gff" -g "$fasta" | \
    awk '$3 ~ /mRNA|gene/ {print $1,$3,$4,$5,$9}'
}

while read line
    do
    fields=($(echo $line | awk '! /#/ && $2 ~ /Complete/ {print $1,$3,$4,$5,$6,$10}'))
        case ${fields[4]} in
        '-')
            region=${fields[1]}:${fields[3]}..${fields[2]};;
        '+')
            region=${fields[1]}:${fields[2]}..${fields[3]};;
        *)
            continue;;
    esac
    on_gff=($(read_region))
    ids="${fields[0]}\t${on_gff[1]}\t${on_gff[0]}\t"
    locs="${on_gff[2]}\t${on_gff[3]}\t"
    info="${on_gff[4]}\t${fields[5]}"
    echo -e "${ids}${locs}${info}" >> gff_busco.txt
    done < "$busco_ref"


# fields[10] busco descriptor
#
#   read_region awk fields
# $1 contig id
# $3 the type of annotation (gene or mRNA)
# $4 the start of the annotation
# $5 the end
# $9 ncbi descriptor

# on_gff[0] contig id
# on_gff[1] annotation type (gene or mRNA)
# on_gff[2] annotation start
# on_gff[3] annotation stop
# on_gff[4] ncbi descriptor
