#!/bin/bash
#SBATCH -J DEPICT
#SBATCH -n NUMBER_OF_THREADS
#SBATCH -o depict.out
#SBATCH -e depict.err
#SBATCH -t 3-0:0
#SBATCH --export ALL

$DEPICT/src/python/depict.py depict.cfg
