#!/bin/bash
source ~/.bash_profile

module unload gcc
module load gcc/9.2.0

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate casereport

Rscript $@
