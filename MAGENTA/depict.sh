#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/depict.out
#$ -e $HOME/depict.err
#$ -N depict2
#$ -pe make -5
#$ -q all.q
#$ -t 1

cd /genetics/bin/MAGENTA_software_package_vs2_July2011
# -singleCompThread
matlab -nodisplay -r "depict, exit"
