#!/bin/bash

LC_ALL='en_US.UTF-8'


KER_FNAME='my_ker.txt'
KER_TYPE='0'
KER_SIZE='31'
KER_WGHT='0.5'

MAKEFLAG=0
RUNFLAG=0
MPIFLAG=0

#VALGRINDCMD='valgrind --leak-check=full --show-leak-kinds=all --suppressions=/usr/share/openmpi/openmpi-valgrind.supp'
VALGRINDCMD='valgrind'
PERFCMD='perf stat -e task-clock,cycles,instructions,cache-references,cache-misses,branches,branch-misses'
#PERFCMD='perf record -e task-clock,cycles,instructions,cache-references,cache-misses,branches,branch-misses'
MPICMD='mpirun -np'
MPI_PROCS=0
MAKERULE=''
CFLAGS=''
CMD=''
EXE='./blur.x' 
MPI_EXE='./blur_mpi.x' 
OMP_EXE='./blur_omp.x' 

LOGFLAG=0
LOGFILE="log.txt"

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
		MAKEFLAG=1
		shift
		MAKERULE="$1"
	;;
	-perf|perf)
		CMD="$PERFCMD $CMD"
	;;
	-valgrind|valgrind)
		CMD="$VALGRINDCMD $CMD"
	;;
	-info|info)
		CFLAGS="$CFLAGS -DINFO"
	;;
	-time|time)
		CFLAGS="$CFLAGS -DTIME"
	;;
	-mpi|mpi)
		RUNFLAG=1
		MPIFLAG=1
		shift
		MPI_PROCS=$1
		CMD="$MPI_EXE"
	;;
	-omp|omp)
		RUNFLAG=1
		shift
		export OMP_DYNAMIC=false
		#export OMP_PROC_BIND=cores
		#export OMP_PROC_BIND=spread
		export OMP_WAIT_POLICY=active 
		export OMP_NUM_THREADS=$1
		CMD="$OMP_EXE"
	;;
	-serial|serial)
		RUNFLAG=1
		CMD="$EXE"
	;;
	-log|log)
		LOGFLAG=1
	;;
esac
shift
done
fi

if  [[ $MAKEFLAG -eq 1 ]]
then
	make clean
	MAKERULE="make $MAKERULE CFLAGS='$CFLAGS'"
	printf '\n'
	echo ${MAKERULE}
	printf '\n'
	eval ${MAKERULE}
fi

if [[ $RUNFLAG -eq 1 ]]
then
CMD="$CMD $ARGS"
if [[ $MPIFLAG -eq 1 ]]
then
	CMD="$MPICMD $MPI_PROCS $CMD"
	
fi
fi
if  [[ $LOGFLAG -eq 1 ]]
then
	CMD="$CMD > $LOGFILE"
fi
printf '\n'
echo ${CMD}
printf '\n'
eval ${CMD}

printf '\n'
echo "done"

