#!/bin/bash
#PBS -l nodes=1:ppn=36
#PBS -l walltime=02:00:00
#PBS -q dssc


cd $PBS_O_WORKDIR 
module load   openmpi/4.0.3/gnu/9.3.0

MOVES="100000000000" 



for i in 1 ; do 

DIRNAME=run${i}

mkdir -p  ${DIRNAME}

cp  mpi_pi.x ./${DIRNAME}


cd ${DIRNAME}



for procs in  36 ; do
echo "executing on ", ${procs}, "  processors" 
N_MOVES=$(($MOVES*$procs))
/usr/bin/time mpirun  --mca btl '^openib' -np ${procs} mpi_pi.x  ${N_MOVES} >out.${procs}
done 

cd ..

done
