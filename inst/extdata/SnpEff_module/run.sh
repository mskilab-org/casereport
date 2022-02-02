#!/bin/sh
{
    . ~/.bash_profile

    module load vcftools
    module load bcftools
    module load tabix/1.1
    module unload python
    module load python/2.7.8

    LIBDIR=$1
    REF=$2
    VCF=$3
    OUTDIR=$4
    CONFIG=$5

    echo "Input VCF: $3"

    rm -f input.vcf

    echo "Output directory: $4"

    mkdir -p $OUTDIR

    if [ ${VCF: -3} == ".gz" ]
    then
        echo "gunzipping"
        zcat $VCF > input.vcf
    elif [ ${VCF: -4} == ".bcf" ]; then
        echo "Converting from VCF to BCF"
        bcftools view $VCF > input.vcf &&
        echo "Finished converting from VCF to BCF"
    else
        ln -sf $VCF input.vcf
    fi

    annvcf=$OUTDIR/annotated.vcf
    annbcf=$OUTDIR/annotated.bcf

    #cmd="java -Xmx4g -XX:ParallelGCThreads=1 -jar $LIBDIR/snpEff.jar -c $LIBDIR/snpEff.config  -v $REF input.vcf > $annvcf"
    cmd="java -Xmx4g -XX:ParallelGCThreads=1 -jar $LIBDIR/snpEff.jar -c $CONFIG  -v $REF input.vcf > $annvcf"
    if { [ ! -e $annvcf ] &&
	     echo "$(date) Running $(echo $cmd)" &&
	     eval $cmd; }; then
	echo "completed SnpEff"
    elif [ -e $annvcf ]; then
	echo "SnpEff output VCF present!"
    else
        echo "SnpEff Broke! Removing annotated.vcf, but keeping a copy of it as annotated.vcf.debug for debugging."
        mv annotated.vcf annotated.vcf.debug; exit 1
    fi
    
    #fi

    if [ ! -e $annbcf ]; then
	{ echo "final post-processing..." &&
	      bcftools view $annvcf -O b -o $OUTDIR/annotated.unsorted.bcf && 
	      bcftools sort $OUTDIR/annotated.unsorted.bcf -O b -o $annbcf &&
	      bcftools index $annbcf ; } ||
	    { echo "final postprocessing failed!"; exit 1; }
    else
        echo "Final BCF output already present!"
    fi

    rm -f input.vcf
    rm -f $OUTDIR/annotated.unsorted.bcf
    rm -f annotated.vcf.debug
    echo "done"
}
