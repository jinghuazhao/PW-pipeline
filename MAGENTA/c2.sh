#!/usr/bash

#$ -S /bin/bash
#$ -o /home/jhz22/c2.out
#$ -e /home/jhz22/c2.err
#$ -N MAGENTA_c2
#$ -pe make -5
#$ -q all.q
#$ -l hostname=b01
#$ -t 1

cd /genetics/bin/MAGENTA_software_package_vs2_July2011
# -singleCompThread
matlab -nodisplay -r "c2, exit"
