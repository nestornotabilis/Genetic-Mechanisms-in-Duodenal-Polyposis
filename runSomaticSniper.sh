#!/bin/sh

NORMAL_BAM=$1	### PATH TO NORMAL BAM FILE
TUMOUR_BAM=$2	### PATH TO TUMOUR BAM FILE
REFERENCE=$3	### PATH TO FASTA REFERENCE

# (i) Run somatic sniper
bam-somaticsniper \
        -q 15 \
        -Q 15 \
        -F vcf \
        -f $REFERENCE \
        $TUMOUR_BAM \
        $NORMAL_BAM \
        out.vcf

# (ii) Index
igvtools index out.vcf

