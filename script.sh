#!/bin/sh
#SBATCH --time=1:00
#SBATCH --partition=cpar

module load perf

perf record -g -o my_perf.data ./MD.exe < inputdata1.txt

perf report -i my_perf.data > my_perf_report.txt

perf annotate -i my_perf.data > my_perf_annotate.txt