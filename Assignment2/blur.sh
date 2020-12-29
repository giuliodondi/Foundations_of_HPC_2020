#!/bin/bash

LC_ALL='en_US.UTF-8'
KER_FNAME='my_ker.txt'
KER_TYPE='0'
KER_SIZE='31 3'
KER_WGHT='0.5'

#VALGRINDCMD='valgrind --leak-check=full --show-leak-kinds=all --suppressions=/usr/share/openmpi/openmpi-valgrind.supp'
VALGRINDCMD='valgrind'
PERFCMD='perf stat -e task-clock,cycles,instructions,cache-references,cache-misses'
#PERFCMD='perf record -e cache-misses'
MPICMD='mpirun -np'
MPIFLAG=0
MPI_PROCS=0
CMD=''
EXE='./blur.x' 
MPI_EXE='./blur_mpi.x' 
OMP_EXE='./blur_omp.x' 

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
ARGS="-input ${FNAME} -kernel-type ${KER_TYPE} -kernel-size ${KER_SIZE} -kernel-weight ${KER_WGHT}"
#ARGS="-input ${FNAME} -kernel-file ${KER_FNAME}"
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
		MPIFLAG=1
		shift
		MPI_PROCS=$1
		CMD="$MPI_EXE $ARGS"
	;;
	-omp|omp)
		shift
		export OMP_NUM_THREADS=$1
		export OMP_PROC_BIND=true
		CMD="$OMP_EXE $ARGS"
	;;
	-serial|serial)
		CMD="$EXE $ARGS"
	;;
esac
shift
done
fi
if [[ $MPIFLAG -eq 1 ]]
then
	CMD="$MPICMD $MPI_PROCS $CMD"
fi
echo ${CMD}
eval ${CMD}
