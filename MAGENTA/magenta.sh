#!/bin/bash

#$ -S /bin/bash
#$ -o $HOME/magenta.out
#$ -e $HOME/magenta.err
#$ -pe make -5
#$ -q all.q
#$ -t 1

matlab -nodisplay -r "magenta, exit"
