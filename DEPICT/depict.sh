#!/bin/bash

#$ -S /bin/bash
#$ -o depict.out
#$ -e depict.err
#$ -pe make -NUMBER_OF_THREADS
#$ -q all.q
#$ -t 1

$DEPICT/src/python/depict.py depict.cfg
