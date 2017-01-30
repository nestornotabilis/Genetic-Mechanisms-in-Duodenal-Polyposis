#!/bin/sh

NORMAL_BAM=$1	### PATH TO NORMAL BAM FILE
TUMOUR_BAM=$2	### PATH TO TUMOUR BAM FILE
REFERENCE=$3	### PATH TO FASTA REFERENCE

# (i) Create pileups
samtools view -b -u -q 15 $NORMAL_BAM | \
	samtools mpileup -f $REFERENCE -Q 10 - \
	> normal.mpileup 

samtools view -b -u -q 15 $TUMOUR_BAM | \
	samtools mpileup -f $REFERENCE -Q 10 - > \
	tumour.mpileup

# (ii) Run VarScan
java -Xmx20g -jar VarScan.jar somatic \
	normal.mpileup \
	tumour.mpileup out \
	—output-vcf 1 \
	—strand-filter 1

