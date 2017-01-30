#!/bin/sh 

NORMAL_BAM=$1		### PATH TO NORMAL BAM FILE
TUMOUR=BAM=$2		### PATH TO TUMOUR BAM FILe
REFERENCE=$3 		### PATH TO FASTA REFERENCE
TARGET_REGIONS=$4	### PATH TO BED FILE ###
DBSNP=$5		### PATH TO DBSNP VCF FILE]
PATH_TO_GATK=$6		### PATH TO GATK DISTRIBUTION


# (i) Run UnifiedGenotyper 
java -Xmx40g -jar $PATH_TO_GATK/GenomeAnalysisTK.jar \
        -R $REFERENCE \
        -T UnifiedGenotyper \
        -I $NORMAL \
        -I $TUMOUR \
        --dbsnp $DBSNP \
        -o prefiltered.out.vcf \
        -glm SNP \
        -stand_call_conf 30.0 \
        -stand_emit_conf 10.0 \
        -dcov 10000 \
        -metrics UnifiedGenotyper.metrics \
        --annotation FisherStrand \
        --annotation Coverage \
        --annotation AlleleBalance \
        --annotation ReadPosRankSumTest \
        --annotation MappingQualityRankSumTest \
        --annotation RMSMappingQuality \
        --annotation QualByDepth \
        --num_threads 8 \
        -L $TARGET_REGIONS


# (ii) Run VariantFiltration
java -Xmx40g -jar $PATH_TO_GATK/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R $REFERENCE \
        --variant prefiltered.out.vcf \
        -o out.vcf \
        --clusterWindowSize 10 \
        --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' \
        --filterName 'HARD_TO_VALIDATE' \
        --filterExpression 'DP < 5' \
        --filterName 'LowCoverage' \
        --filterExpression 'QUAL < 30.0' \
        --filterName 'VeryLowQual' \
        --filterExpression 'QUAL > 30.0 && QUAL < 50.0' \
        --filterName 'LowQual' \
        --filterExpression 'QD < 1.5' \
        --filterName 'LowQD' \
        --filterExpression "FS > 100.0 " \
        --filterName "StrandBias" \
        -L $TARGET_REGIONS



