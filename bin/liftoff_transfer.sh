#!/bin/bash
# Transfer annotations from busco_to_gff file to new locations in a lifted-over
#   gff

liftover=$1
busco_gff=$2
output=$3

info="#BUSCO ID\tType\tSequence\tStart\tStop\tNCBI descriptor\tBUSCO descriptor"
echo -e "$info" > "$output"

function match {
    awk -v geneID="$ID" '$9 ~ geneID {print; exit}' "$liftover"
}

while IFS=$'\t' read -r id type seq start stop ncbi_des busco_des
    do
        if [ -z "$seq" ]; then
            continue
        fi
        ID=$(echo "$ncbi_des" | sed 's/.*geneID/ID/')
        IFS=$'\t' read -r contig _ _ l_start l_stop _ <<< "$(match)"
        ids="$id\t$type\t"
        locs="${contig}\t${l_start}\t${l_stop}\t"
        info="${ncbi_des}\t${busco_des}"
        echo -e "${ids}${locs}${info}" >> "$output"
    done < "$busco_gff"

# fields[0] busco id
# fields[1] type
# fields[2] ncbi descriptor
# fields[3] busco descriptor
