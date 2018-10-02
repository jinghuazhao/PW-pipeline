#!/bin/bash

#SBATCH -n NUMBER_OF_THREADS
#SBATCH -o magma.out
#SBATCH -e magma.err
#SBATCH -export ALL

$DEPICT/src/python/depict.py depict.cfg
