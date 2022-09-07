#!/bin/bash

# Memory request for 8G
#$ -l h_vmem=8G

# Cores
#$ -pe smp 4
#$ -binding linear:4

# Single output files (merge std out with std error)
#$ -j y

# Runtime request
#$ -l h_rt=300:00:00

# Output path for std out
#$ -o /xchip/beroukhimlab/siyun/germline_classifier/outputs/bedpe_combination

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
reuse -q R-4.0

##################
### Run script ###
##################

Rscript /xchip/beroukhimlab/siyun/germline_classifier/scripts/bedpe_combination/20220905_vcf2bedpe_andcombine.R /xchip/beroukhimlab/siyun/testing/all_samples.txt
