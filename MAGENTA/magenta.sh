#!/usr/bash

#$ -S /bin/bash
#$ -o $HOME/magenta.out
#$ -e $HOME/magebta.err
#$ -N MAGENTA
#$ -pe make -5
#$ -q all.q
#$ -t 1

cd /genetics/bin/MAGENTA_software_package_vs2_July2011
# -singleCompThread
matlab -nodisplay -r "magenta, exit"
