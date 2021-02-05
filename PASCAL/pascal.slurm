#!/usr/bin/bash

#SBATCH --job-name=PASCAL
#SBATCH --account ACCOUNT
#SBATCH --partition PARTITION
#SBATCH --qos=PARTITION
#SBATCH --mem=28800
#SBATCH --ntasks=NUMBER_OF_THREADS
#SBATCH --time=3-00:00:00
#SBATCH --output=pascal.out
#SBATCH --error=pascal.err
#SBATCH --export ALL

$PASCAL/Pascal --set=settings.txt --pval=vegas2v2 --runpathway=on --mafcutoff=0.05
