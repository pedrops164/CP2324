#!/bin/sh
#SBATCH --time=3:00
#SBATCH --partition=cpar
#SBATCH --cpus-per-task=40


module load perf

perf record -g -o my_perf_seq.data make runseq
perf report -i my_perf_seq.data > my_perf_seq_report.txt
perf annotate -i my_perf_seq.data > my_perf_seq_annotate.txt

perf record -g -o my_perf_par.data make runpar
perf report -i my_perf_par.data > my_perf_par_report.txt
perf annotate -i my_perf_par.data > my_perf_par_annotate.txt