#!/usr/bin/bash

#SBATCH --job-name=DEPICT
#SBATCH --account ACCOUNT
#SBATCH --partition PARTITION
#SBATCH --qos=PARTITION
#SBATCH --mem=28800
#SBATCH --ntasks=NUMBER_OF_THREADS
#SBATCH --time=1-00:00:00
#SBATCH --output=depict.out
#SBATCH --error=depict.err
#SBATCH --export ALL

$DEPICT/src/python/depict.py depict.cfg
