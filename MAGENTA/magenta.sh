#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/magenta.out
#$ -e $HOME/magebta.err
#$ -N MAGENTA
#$ -pe make -5
#$ -q all.q
#$ -t 1

cd $MAGENTA
matlab -nodisplay -r "magenta, exit"
