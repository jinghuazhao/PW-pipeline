#!/bin/bash

#$ -S /bin/bash
#$ -o magma.out
#$ -e magma.err
#$ -pe make -5
#$ -q all.q
#$ -t 1

magma --gene-results ${_db}.genes.raw --set-annot $db self-contained --model fwer --out $db
