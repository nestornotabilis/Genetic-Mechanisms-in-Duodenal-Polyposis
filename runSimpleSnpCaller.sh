#!/bin/sh

NORMAL_BAM=$1	### PATH TO NORMAL BAM FILE
TUMOUR_BAM=$2	### PATH TO TUMOUR BAM FILE
REFERENCE=$3	### PATH TO FASTA REFERENCE

# (i) Create raw mpileup 
samtools mpileup -ABd100000 -f $REFERENCE $NORMAL_BAM > normal.mpileup
samtools mpileup -ABd100000 -f $REFERENCE $TUMOUR_BAM > tumour.mpileup

# (ii) Run script
simpleSnpCaller.pl -i normal.mpileup -s normal -m 2 > normal.vcf 
simpleSnpCaller.pl -i tumour.mpileup -s tumour -m 2 > tumour.vcf 

# (iii) Index resultant vcf files
igvtools index normal.vcf
igvtools index tumour.vcf

# (iv) Merging into a single vcf file.
bgzip -c tumour.vcf > tumour.vcf.gz
bgzip -c normal.vcf > normal.vcf.gz
tabix -p vcf normal.vcf.gz
tabix -p vcf tumour.vcf.gz
vcf-merge normal.vcf.gz tumour.vcf.gz > out.vcf

# (v) Create index for final merged vcf file.
igvtools index out.vcf

