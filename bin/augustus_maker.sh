#!/usr/bin/bash

assembly=$1
gff=$2
name=$3
initialGB=$4
evaluate=$5
species=$6

if [[ -d ${AUGUSTUS_CONFIG_PATH}/species/${species} ]] ; then
    species_exists=true
    else species_exists=false
fi

if species_exists; then
    new_species.pl --species=$species
fi

# Prepare files
gff3_merge -d $gff
awk '{if ($2=="maker") print }' *all.gff > ${name}_filtered.gff
gff2gbSmallDNA.pl ${name}_filtered.gff ${assembly} 2000 ${initialGB}

# Perform initial training
etraining --species=${species} ${initialGB}
randomSplit.pl ${initialGB} ${params.augsplit}

# Evaluate first model
mv ${initialGB}.test ${evaluate}
augustus --species=${species} ${evaluate} >& AUG_first_evaluation.out
randomSplit.pl ${initialGB} 1000

# Optimize model
optimize_augustus.pl --species=${species} --rounds=3 --onlytrain=${name}.gb.train \
${name}.gb.test >& AUG_log
etraining --species=${species} ${initialGB}

# Evaluate for second time
augustus --species=${species} ${evaluate} >& AUG_second_evalution.out
touch 2-${species}
