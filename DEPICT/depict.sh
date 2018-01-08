#!/bin/bash

#$ -S /bin/bash
#$ -o $PWD/depict.out
#$ -e $PWD/depict.err
#$ -pe make -5
#$ -q all.q
#$ -t 1

$DEPICT/depict.py depict.cfg
