SALSA runs in publicly-available docker containers that contain all the necessary dependencies for code execution.

The first container is build on the broadinstitute/gatk image and will generate phased and annotated single cell allele-specific counts from 10X Genomics single cell gene expression, single cell ATAC or single cell Multiome libraries. 

The following dependencies are included in the p4rkerw/salsa:count_1.0 container:

GATK 4.2.0.0
bwa
STAR
bcftools
hdf5
pysam
shapeit 4.2
WASP 0.3.4

This workflow assumes that you have already aligned your raw data with cellranger, cellranger-atac, or cellranger-arc to generate coordinate-sorted bam files. For additional information, please consult the 10X Genomics website: https://www.10xgenomics.com/



(Follow the steps in order) 
1. Genotype the snRNA or snATAC libraries using GATK (or obtain a vcf from another method)
2. (Optional) - Merge genotypes obtained from matched snRNA-snATAC or multimodal libraries
3. Phase the genotype using shapeit and the 1000G reference    
4. (Optional) - Annotate the vcf with GATK Funcotator to evaluate gnomAD MAF and variant context  
5. Filter the 10X genomics positions sorted bam file by cell barcode using the subsetbam utility  
6. Apply the WASP pipeline to barcode-filtered bam files to perform variant-aware realignment  
7. Obtain pseudobulk, cell-specific, or single cell allele-specific counts with GATK ASEReadCounter  
8. Analyze allele-specific expression across multiple samples with ASEP


