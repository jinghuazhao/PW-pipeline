#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/c2.out
#$ -e $HOME/c2.err
#$ -N c2
#$ -pe make -5
#$ -q all.q
#$ -t 1

cd $MAGENTA
matlab -nodisplay -r "c2, exit"
