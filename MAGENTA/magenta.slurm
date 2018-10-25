#!/bin/bash
#SBATCH -J MAGENTA
#SBATCH -n 5
#SBATCH -o magenta.out
#SBATCH -e magenta.err
#SBATCH -p PARTITION
#SBATCH -t 4-0:0
#SBATCH --export ALL

matlab -nodisplay -nodesktop -nosplash -nojvm -r "Run_MAGENTA_vs2_July_2011, exit"

