#!/bin/bash

#$ -S /bin/bash
#$ -o $HOME/depict.out
#$ -e $HOME/depict.err
#$ -pe make -5
#$ -q all.q
#$ -t 1

$DEPICT/depict.py depict.cfg
