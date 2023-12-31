#!/bin/sh

#export OMP_NUM_THREADS=24

#module load perf
module load gcc/11.2.0

#perf stat -e instructions,cycles,cache-references,cache-misses make runpar > my_perf_par_stat.txt
#perf stat -e instructions,cycles,cache-references,cache-misses mpirun -np 2 make runpar > outputmpi.txt
mpirun -np 25 perf stat -e instructions,cycles,cache-references,cache-misses make runmpi > outputmpi.txt