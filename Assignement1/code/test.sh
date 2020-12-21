#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:00:01
#PBS -q dssc

cd $PBS_O_WORKDIR 

MOVES="100000" 


/usr/bin/time ./pi.x 1000000 > out.pitest

