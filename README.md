# Genetic-Mechanisms-in-Duodenal-Polyposis
A collection of scripts for running variant callers to identify somatic mutations as detailed in forthcoming research paper: Thomas et al., Genetic Mechanisms in Duodenal Polyposis.

Analysis was performed by running five separate variant callers then comparing their outputs.  The five variant callers were:

* GATK UnifiedGenotyper (version 3.3.0)
* SAMTools (version 0.1.19)
* somaticSniper (version 1.0.4)
* varScan2 (version 2.4.0)
* simpleSnpCaller (in-house script)

We provide a shell script for each caller to illustrate the parameter settings used in each case:

## runUnifiedGenotyper.sh 

UnifiedGenotype is one of the variant callers provided by the Genome Analysis Toolkit, [GATK] (https://software.broadinstitute.org/gatk/).

The provided shell script `runUnifiedGenotyper.sh` runs the GATK tools **UnifiedGenotyper** followed by **VariantFiltration**. GATK is written in Java and will require Java runtime to be installed co-installed. 

Usage, e.g., 
`./runUnifiedGenotyper.sh \
	normal.bam \
	tumour.bam \
	ucsc.hg19.fasta \ 
	All_Exon_50mb.bed \
	dbsnp_138.hg19.vcf \
	/share/apps/GATK_3.3.0`

Where input parameters are:
* BAM of normal sample.
* BAM of tumour sample.
* Path to human reference hg19 in fasta format.
* Path to bed file describing target regions.
* Path to vcc file representing dbSNP (version 138 used).
* Path of GATK distribution.

**References:** 
* The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303 [Article](http://genome.cshlp.org/content/20/9/1297)
* A framework for variation discovery and genotyping using next-generation DNA sequencing data DePristo M, Banks E, Poplin R, Garimella K, Maguire J, Hartl C, Philippakis A, del Angel G, Rivas MA, Hanna M, McKenna A, Fennell T, Kernytsky A, Sivachenko A, Cibulskis K, Gabriel S, Altshuler D, Daly M, 2011 NATURE GENETICS 43:491-498 [Article](http://www.nature.com/ng/journal/v43/n5/full/ng.806.html)
* From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M, 2013 CURRENT PROTOCOLS IN BIOINFORMATICS 43:11.10.1-11.10.33 [Article](http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bi1110s43/abstract;jsessionid=DED6F298045EEBB98C8FE6539D74C72E.f04t03)


## runSamtools.sh

[Samtools](http://www.htslib.org) is a toolkit of programs for handling and manipulating SAM and BAM format files including variant calling.  

The provided shell script `runSamtools.sh` runs **samtools mpileup** on the supplied data and then **bcftools view** on the resulting output.  A final step uses **igvtools index**, part of the [IGVTools](https://software.broadinstitute.org/software/igv/igvtools) toolkit, to index the resulting vcf outfile. 

e.g., 
`./runSamtools.sh normal.bam tumour.bam ucsc.hg19.fasta`
Where input parameters are:
* BAM of normal sample.
* BAM of tumour sample.
* Path to human reference hg19 in fasta format.

**References**
* Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9.  [PMID: 19505943](https://www.ncbi.nlm.nih.gov/pubmed/19505943)
* Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627](https://www.ncbi.nlm.nih.gov/pubmed/21903627)
* Li H. Improving SNP discovery by base alignment quality. Bioinformatics. 2011 Apr 15;27(8):1157-8. doi: 10.1093/bioinformatics/btr076. Epub 2011 Feb 13. [PMID: 21320865](https://www.ncbi.nlm.nih.gov/pubmed/21320865)

## runSomaticSniper.sh

[SomaticSniper](http://gmt.genome.wustl.edu/packages/somatic-sniper/) is a dedicated somatic caller of single nucleotide variants. 

The provided script `runSomaticSniper.sh` runs **bam-somaticsniper** followed by **igvtools index**, part of the [IGVTools](https://software.broadinstitute.org/software/igv/igvtools) toolkit, to index the resulting vcf outfile.

e.g., 
`./runSomaticSniper.sh normal.bam tumour.bam ucsc.hg19.fasta`
Where input parameters are:
* BAM of normal sample.
* BAM of tumour sample.
* Path to human reference hg19 in fasta format.

**Reference**
* Larson, D.E., et al., (2011) SomaticSniper: identification of somatic point mutations in whole genome sequencing data. Bioinformatics, 28(3):311-317 [Article](https://academic.oup.com/bioinformatics/article/28/3/311/188933/SomaticSniper-identification-of-somatic-point)

## runVarScan.sh

[VarScan2](http://varscan.sourceforge.net) is a toolkit of programs that includes the somatic caller **somatic**.

The provided shell script `runVarScan.sh` first uses [SAMTools](http://www.htslib.org) to create pileup outputs of each BAM, followed by the **VarScan somatic** command to create an indexed vcc outfile. 

e.g., 
`./runVarScan.sh \
	normal.bam \
	tumour.bam \
	ucsc.hg19.fasta \
	/share/apps/`

Where input parameters are:
* BAM of normal sample.
* BAM of tumour sample.
* Path to human reference hg19 in fasta format.
* Path to VarScan.jar file.

**Reference**
* VarScan 1: Koboldt DC, Chen K, Wylie T, Larson DE, McLellan MD, Mardis ER, Weinstock GM, Wilson RK, & Ding L (2009). VarScan: variant detection in massively parallel sequencing of individual and pooled samples. Bioinformatics (Oxford, England), 25 (17), 2283-5 [PMID: 19542151](https://www.ncbi.nlm.nih.gov/pubmed/19542151)
VarScan 2: Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research [DOI: 10.1101/gr.129684.111](http://genome.cshlp.org/content/22/3/568) 

## runSimpleSnpCaller.sh

The provided **simpleSnpCaller.pl** perl script is a simple in-house program for identifying single nucleotide variants with respect to reference. Output is complete list of all variants without filtering (beyond specifying minimum coverage depth) to give a baseline for comparing true variant callers.

The supplied wrapper shell script `runSimpleSnpCaller.sh` uses  [SAMTools](http://www.htslib.org) to create pileup outputs of each BAM.  It then runs **simpleSnpCaller.pl** on each output, followed by **igvtools index**.  Lastly the resulting vcf files are merged into a single file using **bgzip**, **tabix** part of the [htslib](http://www.htslib.org/download/) distribution, followed by vcf-merge, part of [VCFTools](http://vcftools.sourceforge.net).  Finally the resultant vcc file is indexed with **igvtools index**.

e.g.,
`./runSimpleSnpCaller.sh normal.bam tumour.bam ucsc.hg19.fasta` 
Where input parameters are:
* BAM of normal sample.
* BAM of tumour sample.
* Path to human reference hg19 in fasta format.

**References for third party tools used**

*IGVTools*
* James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 24–26 (2011) [Article] (http://www.nature.com/nbt/journal/v29/n1/abs/nbt.1754.html)
* Helga Thorvaldsdóttir, James T. Robinson, Jill P. Mesirov. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration.  Briefings in Bioinformatics 14, 178-192 (2013). [Article](https://academic.oup.com/bib/article-lookup/doi/10.1093/bib/bbs017)

*Samtools*
* Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9.  [PMID: 19505943](https://www.ncbi.nlm.nih.gov/pubmed/19505943)
* Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627](https://www.ncbi.nlm.nih.gov/pubmed/21903627)
* Li H. Improving SNP discovery by base alignment quality. Bioinformatics. 2011 Apr 15;27(8):1157-8. doi: 10.1093/bioinformatics/btr076. Epub 2011 Feb 13. [PMID: 21320865](https://www.ncbi.nlm.nih.gov/pubmed/21320865)