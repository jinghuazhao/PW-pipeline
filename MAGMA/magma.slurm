#!/usr/bin/bash

#SBATCH --job-name=MAGMA
#SBATCH --account ACCOUNT
#SBATCH --partition PARTITION
#SBATCH --qos=PARTITION
#SBATCH --mem=28800
#SBATCH --ntasks=NUMBER_OF_THREADS
#SBATCH --time=3-00:00:00
#SBATCH --output=magma.out
#SBATCH --error=magma.err
#SBATCH --export ALL

magma --gene-results ${_db}.genes.raw --set-annot $db self-contained --model fwer --out $db
