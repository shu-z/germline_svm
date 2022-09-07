#!/bin/bash

# Memory request for 8G
#$ -l h_vmem=8G

# Cores
#$ -pe smp 4
#$ -binding linear:4

# Single output files (merge std out with std error)
#$ -j y

# Runtime request
#$ -l h_rt=5:30:00

# Output path for std out
#$ -o /xchip/beroukhimlab/siyun/germline_classifier/outputs/fuzzy_matching

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
reuse -q R-4.0

##################
### Run script ###
##################

Rscript /xchip/beroukhimlab/siyun/germline_classifier/scripts/fuzzy_matching/20220830_fuzzymatch.R 
