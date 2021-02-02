#!/bin/bash
#SBATCH -A ACCOUNT
#SBATCH -J MAGMA
#SBATCH -n 5
#SBATCH -o magma.out
#SBATCH -e magma.err
#SBATCH -p PARTITION
#SBATCH -t 3-0:0
#SBATCH -export ALL

magma --gene-results ${_db}.genes.raw --set-annot $db self-contained --model fwer --out $db
