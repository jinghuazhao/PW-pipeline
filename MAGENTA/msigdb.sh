#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/msigdb.out
#$ -e $HOME/msigdb.err
#$ -N msigdb
#$ -pe make -5
#$ -q all.q
#$ -t 1

cd $MAGENTA
matlab -nodisplay -r "msigdb, exit"
