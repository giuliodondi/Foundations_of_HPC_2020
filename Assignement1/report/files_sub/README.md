# Assignment 1

**Section 1**

The theoretical models were implemented as functions in MATLAB R2020, which was used to generate scaling plots and calcvulate the maximum scaling values.

**Section 2**

Bash scripts were written to launch all the required runs and save them in their respective directories. 
All runs were submitted on the dssc queue and requested node 1 qith 48 cores, in order to utiise GPU nodes at all times.
The script loaded the module  `openmpi/4.0.3/gnu/9.3.0` and ran the following bash command :

/usr/bin/time mpirun  --mca btl '^openib' -np ${procs} mpi_pi.x  ${N_MOVES} >out.${procs}

where ${N_MOVES} was given and ${procs} was modified in a loop.
Data was collected from the folders and formatted by means of bash scripts.
LibreOffice Calc (a freeware clone of Microsoft Excel available on Ubuntu) was used to gather all the data into .csv files and average the data, which was then imported into the MATLAB workspace for plotting and fitting with a script.
