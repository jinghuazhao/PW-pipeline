#!/bin/bash
#SBATCH -J PASCAL
#SBATCH -n 5
#SBATCH -o pascal.out
#SBATCH -e pascal.err
#SBATCH -p long
#SBATCH -t 3-0:0
#SBATCH --export ALL

$PASCAL/Pascal --set=settings.txt --pval=vegas2v2 --runpathway=on --mafcutoff=0.05
