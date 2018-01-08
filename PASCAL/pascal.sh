#!/bin/bash

#$ -S /bin/bash
#$ -o $PWD/pascal.out
#$ -e $PWD/pascal.err
#$ -pe make -5
#$ -q all.q
#$ -t 1

$PASCAL/Pascal --set=settings.txt --pval=vegas2v2 --runpathway=on --mafcutoff=0.05
