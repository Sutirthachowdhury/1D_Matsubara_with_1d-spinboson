#!/bin/bash
#SBATCH -p exciton -A exciton
#SBATCH -o rpmd.out
#SBATCH --mem-per-cpu=1GB
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#module load lapack/3.7.0/b1

time ./mat_1d.exe
