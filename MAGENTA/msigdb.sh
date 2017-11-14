#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/msigdb.out
#$ -e $HOME/msigdb.err
#$ -N msigdb
#$ -pe make -5
#$ -q all.q
#$ -t 1

cd /genetics/bin/MAGENTA_software_package_vs2_July2011
# -singleCompThread
matlab -nodisplay -r "msigdb, exit"
