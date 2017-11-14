#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/c2.out
#$ -e $HOME/c2.err
#$ -N c2
#$ -pe make -5
#$ -q all.q
#$ -t 1

cd /genetics/bin/MAGENTA_software_package_vs2_July2011
# -singleCompThread
matlab -nodisplay -r "c2, exit"
