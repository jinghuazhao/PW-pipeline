#!/bin/bash
#SBATCH -J DEPICT
#SBATCH -n NUMBER_OF_THREADS
#SBATCH -o depict.out
#SBATCH -e depict.err
#SBATCH -p PARTITION
#SBATCH -t 15:00
#SBATCH --export ALL

$DEPICT/src/python/depict.py depict.cfg
