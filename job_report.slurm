#!/bin/bash

#SBATCH --nodes=1

echo "# JOB SUCCESS REPORT #"

cat $2/job_IDs_$3.txt | while read line
do
    jn=$(echo $line | awk '{print $1}')
    sn=$(echo $line | awk '{print $2}')
    jid=$(echo $line | awk '{print $3}')
    tnj=$(cat $(find $1 -name "*${jid}*out") | wc -l)
    tnjs=$(cat $(find $1 -name "*${jid}*out") | grep success | wc -l)
    awk -v jn=$jn \
        -v sn=$sn \
        -v tnj=$tnj \
        -v tnjs=$tnjs \
        '{if (tnjs/tnj == 1) print jn" jobs ("sn") -- 100% of jobs successful";
          else if (tnjs/tnj == 0) print jn" jobs ("sn") -- FAILED";
          else print jn" jobs ("sn") -- "1-tnjs/tnj*100"% of jobs FAILED; "tnjs/tnj*100"% of jobs successful"}'
done

if [ $4 == "T" ]
then
    rm -r $2
fi

(exit) && echo success
