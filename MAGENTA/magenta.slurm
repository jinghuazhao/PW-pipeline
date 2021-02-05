#!/usr/bin/bash

#SBATCH --job-name=MAGENTA
#SBATCH --account ACCOUNT
#SBATCH --partition PARTITION
#SBATCH --qos=PARTITION
#SBATCH --mem=28800
#SBATCH --time=5-00:00:00
#SBATCH --output=magenta.out
#SBATCH --error=magenta.err
#SBATCH --export ALL

matlab -nodisplay -nodesktop -nosplash -nojvm -r "Run_MAGENTA_vs2_July_2011, exit"

