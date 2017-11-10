#!/usr/bash

#$ -S /bin/bash
#$ -o /scratch/tempjhz22/matlab.out
#$ -e /scratch/tempjhz22/matlab.err
#$ -N MAGENTA
#$ -pe make -5
#$ -q all.q
#$ -l hostname=b07
#$ -t 1

cd /genetics/bin/MAGENTA_software_package_vs2_July2011
# -singleCompThread
matlab -nodisplay -r "depict, exit"
