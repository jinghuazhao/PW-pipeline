#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/depict.out
#$ -e $HOME/depict.err
#$ -N depict2
#$ -pe make -5
#$ -q all.q
#$ -t 1

cd $MAGENT
matlab -nodisplay -r "depict, exit"
