#!/bin/bash

LC_ALL='en_US.UTF-8'

KER_TYPE=0
KER_SIZE=51
KER_WGHT=0.5

VALGRINDCMD='valgrind'
PERFCMD='perf stat -e task-clock,cycles,instructions,cache-references,cache-misses'
MPICMD='mpirun -np'
EXE='./blur.x' 

#check for missing required arguments
if [[ $# -eq 0 ]]
then
        echo "No arguments."
		echo "Usage : " $0 " filename "
		echo "Optional arguments:"
		echo "-make | make : will run make clean, make and then run the program"
		echo "				make can be followed by switches"
		echo "-perf |perf : will run the program under perf -e"
		echo "-valgr | valgr : will run the program under valgrind"
else
FNAME=$1
CMD="$EXE -input ${FNAME} -kernel-type ${KER_TYPE} -kernel-size ${KER_SIZE} -kernel-weight ${KER_WGHT}"
shift
while [[ $# -gt 0 ]]
do
case $1 in 
	-make|make)
		shift
		make clean
		make "$1"
	;;
	-perf|perf)
		CMD="$PERFCMD $CMD"
	;;
	-valgrind|valgrind)
		CMD="$VALGRINDCMD $CMD"
	;;
	-mpi|mpi)
		shift
		echo $1
		CMD="$MPICMD $1 $CMD"
	;;
	-serial|serial)
	;;
esac
shift
done
fi
echo ${CMD}
eval ${CMD}
