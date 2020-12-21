#!/bin/bash
#PBS -l nodes=1:ppn=48
#PBS -l walltime=08:00:00
#PBS -q dssc


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

MOVES="10000000000" 

N_ITER=3


for (( i=1; i<=${N_ITER}; i++ ))
do

DIRNAME=run${i}

mkdir -p  ${DIRNAME}

cp  mpi_pi.x ./${DIRNAME}


cd ${DIRNAME}



for procs in 1 4 8 12 16 20 24 28 32 36 40 44 48 ; do
echo "executing on ", ${procs}, "  processors" 
N_MOVES=$(($MOVES*$procs))
/usr/bin/time mpirun  --mca btl '^openib' -np ${procs} mpi_pi.x  ${N_MOVES} >out.${procs}
done 

cd ..

done
