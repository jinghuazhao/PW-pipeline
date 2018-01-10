#!/bin/bash

#$ -S /bin/bash
#$ -o magenta.out
#$ -e magenta.err
#$ -pe make -5
#$ -q all.q
#$ -t 1

matlab -nodisplay -nodesktop -nosplash -nojvm -r "Run_MAGENTA_vs2_July_2011, exit"

