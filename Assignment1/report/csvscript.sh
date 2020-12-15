#!/bin/bash

LC_ALL='en_US.UTF-8'


N_RUNS=3


#RUNS=('serial' 'parallel1' 'parallel4' 'parallel8' 'parallel12' 'parallel16' 'parallel20' 'parallel24' 'parallel28' 'parallel32' 'parallel36' 'parallel40' 'parallel44' 'parallel48')
RUNS=('parallel1' 'parallel4' 'parallel8' 'parallel12' 'parallel16' 'parallel20' 'parallel24' 'parallel28' 'parallel32' 'parallel36' 'parallel40' 'parallel44' 'parallel48')
#RUNS=( 'parallel1' 'parallel12' 'parallel24' 'parallel48')
#RUNS={ 'parallel28' }

OUTFILE=weak10.csv

MASTERDIR=weak10_2

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

PATHNAME=/home/giulio/Documents/MHPC/Foundations_of_HPC_2020/Assignement1/code/newruns

PATHNAME=${PATHNAME}'/'${MASTERDIR}'/'

TMPFILE=${SCRIPTPATH}'/tmp'
rm ${TMPFILE}

cd 

for filee in ${PATHNAME}*.sh.e* ; do
    grep -oP '(?<=^).*?(?=\(0avgtext)' $filee > ${TMPFILE}
done



OUTFILE=${SCRIPTPATH}'/'${OUTFILE}
rm ${OUTFILE}


awk -v bash_arr="${RUNS[*]}" 'BEGIN{split(bash_arr,arr);i=1;j=1}
{if (i==1) { print ",,,,,\n,run"j",,,,"}; gsub("(user|system|elapsed|%CPU)",",",$0) ; print arr[i]"," $0 ;i++; if (i>length(arr)) {i=1;j++} } ' ${TMPFILE} >> ${OUTFILE}

rm ${TMPFILE}