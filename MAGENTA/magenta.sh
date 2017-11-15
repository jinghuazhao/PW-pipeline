#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/magenta.out
#$ -e $HOME/magenta.err
#$ -N MAGENTA_${db}
#$ -pe make -5
#$ -q all.q
#$ -t 1

matlab -nodisplay -r "magenta, exit"
