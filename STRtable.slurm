#!/bin/bash

#SBATCH --nodes=1

ref=$1
refdir=$(echo ${ref%/*})
refpre=$(echo ${ref%.*})
ec=$2

module load gatk
module load conda
conda activate ${ec}

if [[ ! -f $refpre.dict ]]
then
    echo "Indexing genome (picard)..."
    picard CreateSequenceDictionary -R $1 -O $refpre.dict
fi

if [[ ! -f $refdir/str_table.tsv ]]
then
    echo "Composing STR table..."
    gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' ComposeSTRTableFile \
	 -R $1 \
	 -O $refdir/str_table.tsv
fi

(exit) && echo success
