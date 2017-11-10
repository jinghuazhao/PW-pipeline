#!/usr/bash

#$ -S /bin/bash
#$ -o /home/jhz22/msigdb.out
#$ -e /home/jhz22/msigdb.err
#$ -N MAGNTA_msigdb
#$ -pe make -5
#$ -q all.q
#$ -l hostname=b02
#$ -t 1

cd /genetics/bin/MAGENTA_software_package_vs2_July2011
# -singleCompThread
matlab -nodisplay -r "msigdb, exit"
