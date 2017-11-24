#!/bin/bash

#$ -S /bin/bash
#$ -o $HOME/magma.out
#$ -e $HOME/magma.err
#$ -pe make -5
#$ -q all.q
#$ -t 1

magma --gene-results ${_db}.genes.raw --set-annot $db self-contained --out $(basename $db)
