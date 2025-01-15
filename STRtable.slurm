#!/bin/bash

#SBATCH --job-name=str_table
#SBATCH --nodes=1

ref=$1
refdir=$(echo ${ref%/*})
refpre=$(echo ${ref%.*})
dc=$(cat $2 | awk '{print $1}')
ec=$(cat $2 | awk '{print $2}')

start=`date +%s`
echo $HOSTNAME

aklog
source ${dc}/conda.sh
conda activate ${ec}

echo "Indexing genome (picard)..."
picard CreateSequenceDictionary -R $1 -O $refpre.dict

echo "Composing STR table..."
if [ $4 ]
then
    $3/gatk ComposeSTRTableFile \
	 -R $1 \
	 -O $refdir/str_table.tsv \
	 -L $4
else
    $3/gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' ComposeSTRTableFile \
	 -R $1 \
	 -O $refdir/str_table.tsv
fi

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
