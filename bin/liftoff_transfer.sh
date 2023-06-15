#!/bin/bash
# Transfer annotations from busco_to_gff file to new locations in a lifted-over
#   gff

liftover=$1
busco_gff=$2
output=$3

info="#BUSCO ID\tType\tSequence\tStart\tStop\tNCBI descriptor\tBUSCO descriptor"
echo -e "$info" > "$output"

function match {
    awk -v geneID="$ID" '$9 ~ geneID {print}' "$liftover"
}

while IFS= read -r line
    do
        fields=($(echo $line | awk '! /#/ {print $1,$2,$6,$7}'))
        if [ -z "${fields[2]}" ]; then
            continue
        fi
        ID=$(echo "${fields[2]}" | sed 's/.*geneID/ID/')
        corres=($(match))
        ids="${fields[0]}\t${fields[1]}\t"
        locs="${corres[0]}\t${corres[3]}\t${corres[4]}\t"
        info="${fields[2]}\t${fields[3]}"
        echo -e "${ids}${locs}${info}" >> "$output"
    done < "$busco_gff"

# fields[0] busco id
# fields[1] type
# fields[2] ncbi descriptor
# fields[3] busco descriptor
