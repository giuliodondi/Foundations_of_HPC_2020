#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -q dssc


DIRNAME=weak_images

OUTNAME=weak_
EXT=.pgm

TYPE=0
MAXVAL=65535
BASEXSIZE=8192
BASEYSIZE=8192


mkdir ./${DIRNAME}

for n in 1 4 8 12 16 20 24 28 32 36 40 44 48; do
	
	./generate_pgm.x ${TYPE}

done

