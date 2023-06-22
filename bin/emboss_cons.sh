#!/bin/bash
file=$1

header=$(awk '$0 ~ ">" {print;exit;}' "$file" | sed 's/>S[[:digit:]]\{2\}|/cons|/g')
echo "$header"
name=$(echo "$file" | sed 's/_.*//g')
echo "$name"
cons -sequence "$file" \
    -outseq "${name}"_cons.fasta \
    -name "$header"
