#!/bin/bash

LC_ALL='en_US.UTF-8'

MAKEFLAG=0
MPIFLAG=0
CHKFLAG=0
LOGFLAG=0
LOGFILE="log.txt"


if [ $1 == "-make" ] || [ $1 == "make" ]
then
	MAKEFLAG=1
else
	INNAME="$1"
	OUTNAME="output.pgm"
	CHKNAME="output2.pgm"
fi



KER_FNAME='-kernel-type my_ker.txt'
KER_TYPE='-kernel-type 0'
KER_SIZE='-kernel-size 31'
KER_WGHT='-kernel-weight 0.2'

	

if [[ -n "${INNAME}" ]]
then
	ARGS=" ${ARGS} -input ${INNAME}"
fi
if [[ -n "${OUTNAME}" ]]
then
	ARGS="${ARGS} -output ${OUTNAME}"
fi


ARGS=" ${ARGS} ${KER_TYPE} ${KER_SIZE} ${KER_WGHT}"
#ARGS="-input ${FNAME} -kernel-file ${KER_FNAME}"




#VALGRINDCMD='valgrind --leak-check=full --show-leak-kinds=all --suppressions=/usr/share/openmpi/openmpi-valgrind.supp'
VALGRINDCMD='valgrind'
PERFCMD='perf stat -e task-clock,cycles,instructions,cache-references,cache-misses,branches,branch-misses'
#PERFCMD='perf record -e task-clock,cycles,instructions,cache-references,cache-misses,branches,branch-misses'
MPICMD='mpirun --mca btl "^openib" --map-by core --report-bindings -np'
MPI_PROCS=0
MAKERULE=''
CFLAGS=''
CMD=''
SERIAL_EXE='./blur_serial.x' 
MPI_EXE='./blur_mpi.x' 
OMP_EXE='./blur_omp.x' 


#check for missing required arguments
if [[ $# -eq 0 ]]
then
        echo "No arguments."
		echo "Usage : " $0 " filename "
		echo "Optional arguments:"
		echo "-make | make : will run make with all the specified targets"
		echo "				time and info will trigger compilation flags"
		echo "-perf |perf : will run the program under perf -e"
		echo "-valgr | valgr : will run the program under valgrind"
else

shift

	if  [[ $MAKEFLAG -eq 1 ]]
	then	
		MAKEFLAG=1
		while [[ $# -gt 0 ]]
		do
			case $1 in 
				serial)
					MAKERULE="$MAKERULE $1"
				;;
				serialtest)
					MAKERULE="$MAKERULE $1"
				;;
				omp)
					MAKERULE="$MAKERULE $1"
				;;
				mpi)
					MAKERULE="$MAKERULE $1"
				;;
				check)
					MAKERULE="$MAKERULE $1"
				;;
				check)
					MAKERULE="$MAKERULE $1"
				;;
				-info|info)
					CFLAGS="$CFLAGS -DINFO"
				;;
				-time|time)
					CFLAGS="$CFLAGS -DTIME"
				;;

			esac
		shift
		done
		make clean
		MAKERULE="make $MAKERULE CFLAGS='$CFLAGS'"
		printf '\n'
		echo ${MAKERULE}
		printf '\n'
		eval ${MAKERULE}
else 

	while [[ $# -gt 0 ]]
	do
		case $1 in 
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
			CMD="$MPI_EXE"
		;;
		-omp|omp)
			shift
			export OMP_DYNAMIC=false
			export OMP_PLACES=cores
			#export OMP_PROC_BIND=close
			export OMP_WAIT_POLICY=active 
			export OMP_NUM_THREADS=$1
			CMD="$OMP_EXE"
		;;
		-serial|serial)
			CMD="$SERIAL_EXE"
		;;
		-log|log)
			LOGFLAG=1
		;;
		-check|check)
			CHKFLAG=1
		;;
	esac
	shift
	done
	
	CMD="$CMD $ARGS"
	if [[ $MPIFLAG -eq 1 ]]
	then
		CMD="$MPICMD $MPI_PROCS $CMD"

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

	if [[ CHKFLAG -eq 1 ]]
	then
		make check
		printf '\n'
		echo "Evaluating ${OUTNAME} against ${CHKNAME}"
		CMD="./check_pgm.x ${OUTNAME} ${CHKNAME}"
		eval ${CMD}
		printf '\n'

	fi
	echo "done"
	
	
fi
fi





