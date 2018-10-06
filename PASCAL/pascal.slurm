#!/bin/bash

#SBATCH -n 5
#SBATCH -o pascal.out
#SBATCH -e pascal.err
#SBATCH --export ALL

$PASCAL/Pascal --set=settings.txt --pval=vegas2v2 --runpathway=on --mafcutoff=0.05
