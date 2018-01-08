#!/bin/bash

#$ -S /bin/bash
#$ -o $PWD/magenta.out
#$ -e $PWD/magenta.err
#$ -pe make -5
#$ -q all.q
#$ -t 1

matlab -nodisplay -nodesktop -nosplash -nojvm -r "Run_MAGENTA_vs2_June_2011, exit"

