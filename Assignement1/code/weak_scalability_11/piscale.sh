#!/bin/bash
#PBS -l nodes=1:ppn=48
#PBS -l walltime=10:00:00
#PBS -q dssc


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

MOVES="100000000000" 

N_ITER=3


for (( i=1; i<=${N_ITER}; i++ ))
do

DIRNAME=run${i}

mkdir -p  ${DIRNAME}

cp pi.x mpi_pi.x ./${DIRNAME}


cd ${DIRNAME}


for procs in 1 12 24 48 ; do
echo "executing on ", ${procs}, "  processors" 
P_MOVES=$(($procs*$MOVES))
mpirun  --mca btl '^openib' -np ${procs} mpi_pi.x  ${P_MOVES} >out.${procs}
done 

cd ..

done
