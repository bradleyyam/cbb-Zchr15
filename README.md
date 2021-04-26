# cbb-Zchr15

[Goal] 
We want to find the top 10 genes in Z’s 15th chromosome that carries the highest mutational burden, defined by the assignment as the genes with the highest number of point mutations. 

[Process]
Steps:
The Z SNP VCF file contains all of the point mutations across all of Z’s chromosomes. They have additional information on alleles and genotype data, as well as quality of reads and other information. 

We produce a new VCF file that only contains SNPs from Z’s 15th chromosome that are of acceptable read quality with the snippet of R code below.


We uploaded this new VCF file to VEP which helped us to annotate each of the SNPs in the VCF file by translating their coordinate position to specific gene designations. This annotated VCF file is available for download as a VCF file or as a tab-delimited table. We opted for the tab-delimited table because it was easier to read.
Variant Effect Predictor - Homo_sapiens - GRCh37 Archive browser 103 (ensembl.org)
Click “New Job”
Ensure it is for GRCh37.p13 reference genome
Upload filtered vcf file
Use Ensembl/GENCODE transcripts
Use otherwise default settings
After job has run, select download as “TXT”


We ingested the tab-delimited table of each observation of SNP across many genes and simply tabulated the number of counts each gene was recorded with an SNP. This final count accords with the simple “mutational burden” heuristic of genes with the highest number of point mutations. We take the top 10 identifiable genes for further investigation.
