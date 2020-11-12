#!/bin/bash

LC_ALL='en_US.UTF-8'



RUNS=(1 4 8 12 16 20 24 28 32 36 40 44 48)
#RUNS=(40 44)
#RUNS=(1 12 24 36 48)

OUTFILE=walltimes10.csv


MASTERDIR=weak10_2

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

PATHNAME=/home/giulio/Documents/MHPC/Foundations_of_HPC_2020/Assignement1/code/newruns

PATHNAME=${PATHNAME}'/'${MASTERDIR}'/'

OUTFILE=${SCRIPTPATH}'/'${OUTFILE}
rm ${OUTFILE}



for run in  1 2 3 
do 

	DIR=${PATHNAME}'run'${run}
	
	cd ${DIR}
	
	
	echo ",,,,," >> ${OUTFILE}
	echo ",run"${run}",,,," >> ${OUTFILE}


	for p in ${RUNS[@]}
	do
		echo "run "$run $p
		INFILE=${DIR}'/out.'${p}
		
		#awk -v bash_arr="${RUNS[*]}" 'BEGIN{split(bash_arr,arr);i=1;j=1}
		#{if (i==1) { print ",,,,,\n,run"j",,,,"}; gsub("(user|system|elapsed|%CPU)",",",$0) ; print arr[i]"," $0 ;i++; if (i>length(arr)) {i=1;j++} } ' ${TMPFILE} >> ${OUTFILE}

		##awk -v ' BEGIN{match=0;average=0;maxm=0;minm=0}
		##/(walltime*)/  { match = `grep -oP '[0-9]\+' $0 `; average=(average + match); if (maxm==0 || match>maxm) {maxm=match}; if (minm==0 || match<minm) {minm=match} }
		##END { AVG=(average/$p); MIN=minm;MAX=maxm }' ${INFILE}
		
		# g
		#"grep -oP '(?=:).*'
		#
		
		#| grep -o '0\.[0-9]\+'
		#echo ` awk '/walltime*/ {print $0}' ${INFILE} | grep '[0-9]\+' `
		
		# | grep -o '[0-9\.]\+'
		
		echo -n ${p}','>> ${OUTFILE}
		for line in  `grep walltime ${INFILE}   `
		do
			
			NUM=$(echo $line | grep -o '[0-9]\+\.[0-9]\+' )
			if [ ! -z $NUM ]
			then
				echo -n $NUM',' >> ${OUTFILE}
			fi
			
		done
		echo " " >> ${OUTFILE}
	done
	
done