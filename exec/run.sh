#!/bin/bash
{
    . ~/.bash_profile
    set -a
    
    module unload gcc
    module load gcc/9.2.0

    module load anaconda3
    . /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
    # conda activate casereport 

    set -x
    export condaenv=${1}; shift
    ## "/gpfs/commons/groups/imielinski_lab/Software/anaconda3/envs/casereport"
    set +x

    if [ -e ${condaenv} ] && [ ! "${condaenv}" = "/dev/null" ]; then
	echo "Activated conda environment: $condaenv"
	conda activate ${condaenv}
    else
	echo "$condaenv not a proper conda environment path"
	conda activate casereport
    fi
    
    module load bcftools
    
    cmd="Rscript $@"
    echo "Running R script:" && echo "$cmd" && eval $cmd
    cmdsig=$?

    if [ "$cmdsig" = "0" ]; then
	echo "Completed WGS Case Report!"
    else
	echo "WGS Case Report job failed!"
	exit $cmdsig
    fi

    exit 0

}
