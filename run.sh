#!/bin/bash
source ~/.bash_profile

module unload gcc
module load gcc/9.2.0

Rscript $@
