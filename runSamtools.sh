#!/bin/sh

NORMAL_BAM=$1	### PATH TO NORMAL BAM FILE
TUMOUR_BAM=$2	### PATH TO TUMOUR BAM FILE
REFERENCE=$3	### PATH TO FASTA REFERENCE


# (i) Run pileup
samtools mpileup \
        -DSu \
        -q 15 \
        -Q 10 \
        -C 50 \
        -m 3 \
        -F 0.0002 \
        -d 100000 \
        -f $REFERENCE \
        $NORMAL_BAM \
        $TUMOUR_BAM \
        > out.mpileup

# (ii) Run bcftools
bcftools view -vcgT pair out.mpileup > out.vcf

# (iii) index
igvtools index out.vcf

