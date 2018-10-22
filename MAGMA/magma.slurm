#!/bin/bash
#SBATCH -J MAGMA
#SBATCH -n 5
#SBATCH -o magma.out
#SBATCH -e magma.err
#SBATCH -t 3-0:0
#SBATCH -export ALL

magma --gene-results ${_db}.genes.raw --set-annot $db self-contained --model fwer --out $db
