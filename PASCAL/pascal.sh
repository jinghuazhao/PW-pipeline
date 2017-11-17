#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/pascal.out
#$ -e $HOME/pascal.err
#$ -pe make -5
#$ -q all.q
#$ -t 1

$PASCAL/Pascal --set=settings.txt --pval=vegas2v2 --runpathway=on --mafcutoff=0.05
