#!/bin/bash

#$ -S /bin/bash
#$ -o depict.out
#$ -e depict.err
#$ -pe make -NUMBER_OF_THREADS
#$ -q all.q
#$ -t 1

$DEPICT/depict.py depict.cfg
