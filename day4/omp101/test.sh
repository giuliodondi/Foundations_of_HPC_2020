#PBS	-l nodes=1:ppn=4,walltime=00:05:00
#PBS	-q dssc

pwd
cd $PBS_O_WORKDIR
pwd

./a.out

export OMP_NUM_THREADS=20

./a.out
exit
