# SNPLookUp

SNPLookUp automates and optimises the process of identifying single nucleotide polymorphisms(SNPs) with high Global Minor Allele Frequency(GMAF) and low linkage disequilibrium(LD) within user-defined regions of interest.

SNPLookUp takes as input a bed file containing regions of interest in the GRCh37 chromosome. For each region of interest, SNPs within the region are identified and filtered based on user-defined cutoffs for minimum GMAF, maximum pairwise LD between SNPs within the region, and/or maximum number of SNPs.

SNPs selected can be subsequently visualised and saved in either csv or bed format that are compatible as inputs for Illumina's DesignStudio web service.

### Setup

Python package dependencies:
-	csv
-	gzip
-	matplotlib
-	requests
-	tqdm

Install via “pip install csv gzip matplotlib requests tqdm” on the command line.

In addition, the program requires data from GRCh37p13 dbSNP chromosome reports in order to work for GRCh37 chromosomal positions. retrievable from:
ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/chr_rpts/
Downloaded chromosome reports should be stored in /SNPLookUp/data/snps_by_chr_hg19/.


### Usage

Details about main functions as well as sample usage for this program can be found in sample_run.ipynb

