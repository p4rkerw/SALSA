🌶️SALSA is a tool for generating and analyzing phased single cell allele-specific read counts from 10X Genomics cellranger datasets. For additional information, please consult the 10X Genomics: https://www.10xgenomics.com/ . The workflow runs in a publicly-available docker container with all the necessary dependencies for code execution. For this tutorial, we will download a dataset from the 10X Genomics website that has already been aligned and annotated. To complete all of the steps in the tutorial, you will need to download a cellranger reference, GATK bundle resources, and 1000G phased reference (see below for more information). However, all of the tutorial outputs are included in the this repository in case you want to compare results or skip a step. If you are analyzing your own data you will need to run cellranger and annotate your barcodes before starting the workflow.

**Step 0: Download cellranger reference** If you don't already have a GRCh38 cellranger reference download one from the [10X Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). 10X Genomics routinely updates their references with each new cellranger build, but new references are often backwards-compatible. The cellranger reference in this tutorial is compatible with the tutorial dataset (which is aligned to GRCh38-2020-A) Feel free to download a different reference and/or [dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets) from the 10X Genomics collection; just make sure it's aligned to GRCh38 so it matches the GATK bundle resources and ucsc contig style. Alignment information can be found in the summary html files. 

**Step 0: Download tutorial dataset and align to your chosen reference with cellranger** We will download a a single cell gene expression dataset obtained from 1k PBMCs from a healthy donor and align it to our reference. This dataset uses the single cell gene expression v3 chemistry. 
```
# URL to the dataset: https://support.10xgenomics.com/single-cell-gene-expression/datasets/6.0.0/1k_PBMCs_TotalSeq_B_3p_LT
# create your salsa tutorial directory and download the fastq
project=$PWD/salsa
wget -P $project/tar https://cf.10xgenomics.com/samples/cell-exp/6.0.0/1k_PBMCs_TotalSeq_B_3p_LT/1k_PBMCs_TotalSeq_B_3p_LT_fastqs.tar

# check md5
md5sum $project/tar/1k_PBMCs_TotalSeq_B_3p_LT_fastqs.tar #ac98de1046df421ff3d8dc6b1a3d6112

# unpack
tar -C $project -xvf $project/tar/1k_PBMCs_TotalSeq_B_3p_LT_fastqs.tar

# align and count with cellranger (version 6.0.1)
reference=/g/reference
cellranger count \
--id pbmc_1k \
--fastqs $project/1k_PBMCs_TotalSeq_B_3p_LT_fastqs/1k_PBMCs_TotalSeq_B_3p_LT_gex_fastqs \
--transcriptome $reference/refdata-gex-GRCh38-2020-A \
--nosecondary \
--nopreflight \
--expect-cells 1000 \
--localcores 10
```

**Step 0: Download GATK resource bundle** The following files are required for genotyping with GATK using GRCh38 and can be found in the [GATK google cloud bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false) . Download these files to your reference directory. For additional information on GATK germline and RNA-seq short variant discovery check out their [website](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
```
resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz
resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz
resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz
resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
resources_broad_hg38_v0_wgs_calling_regions.hg38
```
The following additional files are required if you are analyzing single cell ATAC datasets (not needed for the tutorial):
```
resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz
```

**Step 0: Pull 🌶️SALSA container** 
```
docker pull p4rkerw/salsa:latest
```
The 🌶️SALSA Docker container is built on [broadinstitute/gatk:4.2.0.0](https://hub.docker.com/r/broadinstitute/gatk) with additional dependencies pre-installed (see below). Docker is a set of platform as a service (PaaS) products that use OS-level virtualization to deliver software in packages called containers, which can be used to enhance reproducibility. You can find the the 🌶️SALSA Dockerfile [here](https://github.com/p4rkerw/SALSA/blob/main/docker/Dockerfile).
```
# container dependencies
GATK 4.2.0.0
bwa 0.7.17
STAR 2.7.4a
bcftools 1.9 
pysam 0.15.3
shapeit 4.2
WASP 0.3.4
```

**Step 0: Clone 🌶️SALSA github repository** The repository is cloned to the $project directory, which is the same directory that the tutorial dataset was downloaded to. The tutorial assumes that the path to your project directory is /g/salsa so make sure to change the path if you chose a different directory.
```
git -C $project clone https://github.com/p4rkerw/SALSA
```

**Step 1: Genotype a single cell gene expression dataset** The tutorial workflow is based on the GATK germline short variant discovery pipeline for RNAseq. Additional info can be found on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) . To explore additional GATK options type 'gatk --list' into the terminal. If you want to skip ahead move the [genotyped vcf](https://github.com/p4rkerw/SALSA/blob/main/Tutorial/pbmc.rna.chr22.vcf.gz) and its index to the volume mounted to vcfdir/rna_genotype and proceed to the next step.
```
Usage: step1_gatk_genotype.sh [-indomlt]
   -i  | --bam                STR   path/to/input.bam eg. [rna_counts/sample_1/outs/possorted*.bam]
   -n  | --library_id         STR   library_id: eg. [sample_1]
   -d  | --outputdir          STR   output directory name eg. [vcfdir/rna_genotype]
   -o  | --outputvcf          STR   name of output vcf eg. [sample_1.rna.vcf.gz]
   -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
   -l  | --interval           STR   optional: genotype a single chromosome eg. [chr22]
   -t  | --threads            INT   number of threads. Default=[1]
   -h  | --help                     show usage
```
**Launch 🌶️SALSA container** : Mount the required volumes in an interactive session. The $SCRATCH1 variable designates a temporary file directory. The directory with coordinate-sorted bam file generated by cellranger count is mounted to rna_counts. Users may prefer to mount all the directories required for future steps rather than launching a new container for each step (see below), but for the sake of simplicity only the volumes required for each step are specified.
```
SCRATCH1=/g/scratch
project=/g/salsa
reference=/g/reference
docker run \
--workdir $HOME \
-v $HOME:$HOME \
-v $project/pbmc:$HOME/rna_counts \
-v $reference/refdata-gex-GRCh38-2020-A:$HOME/rna_ref \
-v $project/vcf_output:$HOME/vcfdir \
-v $project/SALSA:$HOME/SALSA \
-v $reference/gatk:$HOME/gatk_bundle \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/g/scratch" \
--rm -it p4rkerw/salsa:latest
```
**Genotype an RNA sample** 
```
bash SALSA/step1_gatk_genotype.sh \
--bam rna_counts/outs/possorted_genome_bam.bam \
--library_id pbmc_1k \
--outputdir vcfdir/rna_genotype \
--outputvcf pbmc.rna.chr22.vcf.gz \
--interval chr22 \
--modality rna \
--threads 10
```
**Inspect the genotyped vcf** Note how the GT field has a "/" character, which indicates that the variants are not phased. 
```
bcftools query -f '[%CHROM,%POS,%REF,%ALT,%GT,%FILTER\n]' vcfdir/rna_genotype/pbmc.rna.chr22.vcf.gz | head -n5
# chr22,17085042,T,C,1/1
# chr22,17085084,G,C,1/1
# chr22,17110243,A,G,0/1
# chr22,17110377,A,G,0/1
# chr22,17111083,G,T,0/1
```
**(Not required for tutorial) Step 2: Merge genotypes from the same patient** If you genotyped a paired single cell gene expression and ATAC dataset from a split sample (or a single cell Multiome) you can merge genotypes into a single vcf. If you're following the tutorial, you can skip this step.
```
Usage: step2_merge_geno.sh [-nabdit]
   -n  | --library_id         STR   library_id: eg. [sample_1]
   -a  | --vcfone             STR   path/to/first.vcf.gz eg. [vcfdir/atac_genotype/sample_1.atac.vcf.gz]. Prioritize annotations from first vcf if there is overlap.
   -b  | --vcftwo             STR   path/to/second.vcf.gz eg. [vcfdir/rna_genotype/sample_1.rna.vcf.gz]
   -d  | --outputdir          STR   output directory [vcfdir/joint_genotype]
   -o  | --outputvcf          STR   name of output vcf eg. [sample_1.pass.joint.vcf.gz]
   -i  | --includefiltered          optional: include filtered variants in output vcf. Default=[false]
   -t  | --threads            INT   number of threads. Default=[1]
   -h  | --help                     show usage
```
**Launch 🌶️SALSA container**
```
# SCRATCH1=/g/scratch
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v $project/vcf_output:$HOME/vcfdir \
# -v $project/SALSA:$HOME/SALSA \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/g/scratch" \
# --rm -it p4rkerw/salsa:latest
```
**Merge two genotypes**
```
# bash SALSA/step2_merge_geno.sh \
# --library_id sample_1 \
# --vcfone vcfdir/$atac_genotype/sample_1.atac.chr22.vcf.gz \
# --vcftwo vcfdir/${modalitytwo}_genotype/sample_1.rna.chr22.vcf.gz \
# --outputdir vcfdir/joint_genotype \
# --outputvcf sample_1.pass.joint.chr22.vcf.gz \
# --threads 4
```
**(Recommended) Step 3: Phase genotype** If you want to perform your analysis with phased genotypes you will need a phased reference. We will use [shapeit4](https://github.com/odelaneau/shapeit4) to phase our variants. To explore additional phasing options type 'shapeit4.2' into your terminal. Variant phasing increases the performance of the WASP variant-realignment and downstream analysis steps. Download the 1000G phased reference files for SNV only or SNV_and_INDEL from ftp.1000genomes.ebi.ac.uk . If you are analyzing the tutorial RNA dataset select the SNV reference. If you want to skip ahead move the [phased vcf](https://github.com/p4rkerw/SALSA/blob/main/Tutorial/pbmc.pass.joint.chr22hcphase.vcf.gz) and its index to the volume mounted to vcfdir/phasing and proceed to the next step.

a) SNV only: [/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/) </br>
b) SNV_and_INDEL: [/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/)
```
Usage: step3_phase_vcf.sh [-nvdolpsitrh]
   -n  | --library_id         STR   library_id: eg. [sample_1]
   -v  | --inputvcf           STR   path/to/input.vcf.gz eg. [vcfdir/joint_genotype/sample_1.pass.joint.vcf.gz]
   -d  | --outputdir          STR   output directory name eg. [vcfdir/phasing]
   -o  | --outputvcf          STR   name of output vcf eg. [sample_1.pass.joint.phase.vcf.gz]
   -l  | --interval           STR   optional: phase a single chromosome eg. [chr22]
   -p  | --hcphase                  optional: recover haplotypecaller physical phasing variants that are not in shapeit reference. Default=[false]
   -s  | --snvonly            STR   use the biallelic_SNV reference for phasing
   -i  | --snvindel           STR   use the biallelic_SNV_and_INDEL reference for phasing
   -r  | --reproduce                optional: run shapeit with a single thread for reproducibility. Default=[false]
   -t  | --threads            INT   number of threads. Default=[1]
   -h  | --help                     show usage
```

You will eventually need to download the vcf for every chromosome, but for the purposes of the tutorial just download the SNV reference for chr22 to $reference/phasing/biallelic_SNV . The reference may take awhile to download so feel free to move ahead to the next step using the [phased vcf](https://github.com/p4rkerw/SALSA/blob/main/Tutorial/pbmc.pass.joint.chr22hcphase.vcf.gz) in the repository. You can come back to step 3 when the download finishes.
```
ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
```

**Launch 🌶️SALSA container**
```
SCRATCH1=/g/scratch
project=/g/salsa
reference=/g/reference
docker run \
--workdir $HOME \
-v $HOME:$HOME \
-v $project/vcf_output:$HOME/vcfdir \
-v $project/SALSA:$HOME/SALSA \
-v $reference/phasing/biallelic_SNV:$HOME/phasing \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/g/scratch" \
--rm -it p4rkerw/salsa:latest
```
**Rename the reference contigs** The 1000G vcf reference files do not have the same contig style as the cellranger reference. You will need to update the 1000G contig style using bcftools. For the tutorial, we will only do chromosome 22. 
```
# runtime ~10min
# rename the contigs in the 1000G reference from 22 to chr22
for i in {1..22} X;do echo "${i} chr${i}";done > /tmp/rename_chrm.txt
inputvcf=ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
bcftools annotate phasing/$inputvcf --threads 4 --rename-chrs /tmp/rename_chrm.txt -Oz -o $SCRATCH1/$inputvcf
mv $SCRATCH1/$inputvcf phasing/
bcftools index --threads 4 phasing/$inputvcf
```
**Phase an interval** For the tutorial, we will not set the --reproduce flag so we can use multithreading. As a result, there may be small differences between the phased vcf in the repository and your results. 
```
bash SALSA/step3_phase_vcf.sh \
--library_id pbmc_1k \
--inputvcf vcfdir/rna_genotype/pbmc.rna.chr22.vcf.gz \
--outputdir vcfdir/phasing \
--outputvcf pbmc.pass.joint.chr22hcphase.vcf.gz \
--interval chr22 \
--hcphase \
--snvonly \
--threads 10
```
**Inspect the phased vcf** The GT field has a "|" character, indicating that the variants are now phased.
```
bcftools query -f '[%CHROM,%POS,%REF,%ALT,%GT\n]' vcfdir/phasing/pbmc.pass.joint.chr22hcphase.vcf.gz | head -n5
# chr22,17085042,T,C,1|1
# chr22,17085084,G,C,1|1
# chr22,17111083,G,T,0|1
# chr22,17115288,T,C,0|1
# chr22,17115498,G,C,0|1
```
**(Optional) Step 4: Annotate vcf with GATK Funcotator** If you want to annotate your vcf with GENCODE and gnomAD you will need to download the [GATK Funcotator resource](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial). GATK routinely updates its resources so you may need to change the name of the folder in the tutorial to match the one you downloaded. gnomAD resources need to be enabled after download (see GATK instructions on their website). When the resources have been downloaded move the dataSources folder into to the reference directory (eg. [reference/funcotator_dataSources.v1.6.20190124]). If you want to skip ahead while these files are downloading move the [funcotated vcf](https://github.com/p4rkerw/SALSA/blob/main/Tutorial/pbmc.pass.joint.chr22hcphase.funco.vcf.gz) and its index to the volume mounted to vcfdir/funcotation and proceed to the next step.
```
Usage: step4_gatk_anno_vcf.sh [-nvdoamfth]
   -n  | --library_id         STR   library_id: eg. [sample_1]
   -v  | --inputvcf           STR   path/to/input.vcf.gz eg. [vcfdir/phasing/sample_1.pass.joint.hcphase.vcf.gz]
   -d  | --outputdir          STR   output directory name eg. [vcfdir/funcotation]
   -o  | --outputvcf          STR   name of output vcf eg. [sample_1.pass.joint.hcphase.funco.vcf.gz]
   -a  | --output_table       STR   name of output funcotation csv eg. [sample_1.pass.joint.hchpase.formatted.csv]
   -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
   -f  | --funcotation        STR   path/to/funcotation directory eg. [reference/funcotator_dataSources.v1.6.20190124g]
   -t  | --threads            INT   number of threads. Default=[1]
   -h  | --help                     show usage
```
**Launch 🌶️SALSA container**
```
SCRATCH1=/g/scratch
project=/g/salsa
reference=/g/reference
docker run \
--workdir $HOME \
-v $HOME:$HOME \
-v $reference/refdata-gex-GRCh38-2020-A:$HOME/rna_ref \
-v $project/vcf_output:$HOME/vcfdir \
-v $project/SALSA:$HOME/SALSA \
-v $reference:$HOME/reference \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/g/scratch" \
--rm -it p4rkerw/salsa:latest
```
**Annotate a vcf**
```
# runtime ~3min
bash SALSA/step4_gatk_anno_vcf.sh \
--library_id pbmc_1k \
--inputvcf vcfdir/phasing/pbmc.pass.joint.chr22hcphase.vcf.gz \
--outputdir vcfdir/funcotation \
--outputvcf pbmc.pass.joint.chr22hcphase.funco.vcf.gz \
--output_table pbmc.pass.joint.chr22hcphase.formatted.csv \
--modality rna \
--funcotation reference/funcotator_dataSources.v1.6.20190124g \
--threads 10
```
**Inspect the annotation table** GATK Funcotator provides a lot of annotation fields in the vcf INFO field and only a subset are included in this table. Note how these two variants are present in the gnomAD database and have allele frequency annotations in the gnomAD_exome_AF and gnomAD_genome_AF columns.
```
head -n3 vcfdir/funcotation/pbmc.pass.joint.chr22hcphase.formatted.csv
# variant_id,CHROM,POS,REF,ALT,GT,FILTER,Gencode_27_variantClassification,Gencode_27_codonChange,gnomAD_exome_AF,gnomAD_genome_AF,Gencode_27_hugoSymbol
# chr22_17085042_T_C,chr22,17085042,T,C,1|1,.,FIVE_PRIME_UTR,,8.02559e-01,8.44264e-01,IL17RA
# chr22_17085084_G_C,chr22,17085084,G,C,1|1,.,FIVE_PRIME_UTR,,8.03628e-01,8.44332e-01,IL17RA
```

**(Recommended) Step 5:** Use barcode celltype annotations to filter the coordinate-sorted cellranger bam using the CB tag. This step will speed up downstream analysis by eliminating barcodes that do not meet quality control. The barcode annotation file should have three columns where the first column is the barcode, the second column is the library_id, and the third column is the celltype annotation. For the purposes of the tutorial, we will only filter chr22.
```
Usage: step5_filterbam.sh [-nidolmbeth]
   -n  | --library_id         STR   library_id: eg. [sample_1]
   -i  | --inputbam           STR   path/to/input.bam eg. [rna_counts/sample_1/outs/possorted_genome_bam.bam]
   -d  | --outputdir          STR   output directory eg. [project/wasp_rna]
   -o  | --outputbam          STR   filtered output bam eg. [sample_1.bcfilter.bam]
   -l  | --interval           STR   optional: filter a single chromosome eg. [chr22]
   -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
   -b  | --barcodes           STR   path/to/barcodes.csv with headers and three columns. First column is named "barcodes"
                                    second column is group "orig.ident" eg. [barcodes/rna_barcodes.csv]
   -e  | --validate                 validate the barcode-filtered bam file. Default=[false]
   -t  | --threads            INT   number of threads. Default=[1]
   -h  | --help                     show usage
```
**Launch 🌶️SALSA container**
```
SCRATCH1=/g/scratch
project=/g/salsa
docker run \
--workdir $HOME \
-v $HOME:$HOME \
-v $project/pbmc_1k:$HOME/rna_counts \
-v $project/SALSA:$HOME/SALSA \
-v $project/barcodes:$HOME/barcodes \
-v $project:$HOME/project \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/g/scratch" \
--rm -it p4rkerw/salsa:latest
```
**Download clustering analysis for tutorial dataset and create a barcode csv**
```
# download the barcode cluster annotation file from 10X Genomics
wget -P $project https://cf.10xgenomics.com/samples/cell-exp/6.0.0/1k_PBMCs_TotalSeq_B_3p_LT/1k_PBMCs_TotalSeq_B_3p_LT_analysis.tar.gz
tar -C $project/pbmc_1k -xvzf $project/1k_PBMCs_TotalSeq_B_3p_LT_analysis.tar.gz

# use the cluster number as celltype and assign pbmc as the library_id in final csv
cluster=$project/pbmc_1k/analysis/clustering/graphclust/clusters.csv
(echo "barcode,orig.ident,celltype"; awk 'BEGIN{FS=OFS=","} {if (NR!=1) print $1,"pbmc_1k",$2}' $cluster) |sed 's/-1//g'> barcodes/rna_barcodes.csv
```
**Inspect the barcode csv**
```
head -n3 barcodes/rna_barcodes.csv
# barcode,orig.ident,celltype
# AATCACGAGCAGCCCT,pbmc_1k,1
# AATCACGAGGAACTCG,pbmc_1k,1
```
**Filter cellranger bam with barcode csv**
```
bash SALSA/step5_filterbam.sh \
--library_id pbmc \
--validate \
--inputbam rna_counts/possorted_genome_bam.bam \
--modality rna \
--interval chr22 \
--barcodes barcodes/rna_barcodes.csv \
--outputdir project/wasp_rna \
--outputbam pbmc.bcfilter.chr22.bam \
--threads 10
```
**Step 6: Perform variant-aware realignment with WASP** This step takes a genotyped vcf and performs variant-aware realignment on a coordinate-sorted and indexed bam file with WASP. WASP is a tool to perform unbiased allele-specific read mapping and you can read more about it here: https://github.com/bmvdgeijn/WASP . For the purposes of the tutorial, we will only analyze chromosome 22. For RNA analysis, this step requires a STAR index of the cellranger reference. A STAR index can be built ahead of time using the command below. Building a new index takes awhile, but it only needs to be done once.
```
Usage: step6_wasp.sh [-vbdoginlmpt]
   -v  | --inputvcf          STR   vcfdir/funcotation/sample_1.pass.joint.hcphase.funco.vcf.gz
   -b  | --inputbam          STR   path/to/input.bam eg. [project/wasp_rna/sample_1.bcfilter.bam]
   -d  | --outputdir         STR   name of output directory eg. [project/wasp_rna]
   -o  | --outputbam         STR   name of output wasp bam eg. [sample_1.phase.wasp.bam]
   -g  | --genotype          STR   genotype: [rna] [atac] [joint]
   -i  | --stargenome        STR   path/to/star genomeDir made with genomeGenerate eg. [rna_ref/salsa_star]. If absent create new STAR index in rna_ref/salsa_star
   -n  | --library_id        STR   library_id: eg. [sample_1]
   -m  | --modality          STR   modality: [rna] [atac]
   -l  | --interval          STR   specified interval eg. [chr22]
   -p  | --isphased                input vcf is phased. Default=[FALSE]
   -t  | --threads           INT   number of threads. Default=[1]
   -h  | --help                    show usage
```
**Launch 🌶️SALSA container**
```
SCRATCH1=/g/scratch
project=/g/salsa
reference=/g/reference
docker run \
--workdir $HOME \
-v $HOME:$HOME \
-v $reference/refdata-gex-GRCh38-2020-A:$HOME/rna_ref \
-v $reference:$HOME/reference \
-v $project/vcf_output:$HOME/vcfdir \
-v $project/SALSA:$HOME/SALSA \
-v $project:$HOME/project \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/g/scratch" \
--rm -it p4rkerw/salsa:latest
```

**(Not required for tutorial) Build a STAR index for 🌶️SALSA** The refdata-gex-GRCh38-2020-A reference comes with a pre-packaged STAR reference built with STAR-2.5.1b. If you want to use a more recent version of STAR or a different index you will need to build a new one. 
```
# STAR \
# --runMode genomeGenerate \
# --genomeDir rna_ref/salsa_star \
# --genomeFastaFiles rna_ref/fasta/genome.fa \
# --sjdbGTFfile rna_ref/genes/genes.gtf \
# --genomeSAsparseD 3 \
# --runThreadN 10
```

**Run WASP on the barcode-filtered bam**
```
# runtime ~2min
bash SALSA/step6_wasp.sh \
--inputvcf vcfdir/funcotation/pbmc.pass.joint.chr22hcphase.funco.vcf.gz \
--inputbam project/wasp_rna/pbmc.bcfilter.chr22.bam \
--outputdir project/wasp_rna \
--outputbam pbmc.hcphase.chr22wasp.bam \
--genotype joint \
--stargenome rna_ref/star \
--library_id pbmc_1k \
--modality rna \
--isphased \
--interval chr22 \
--threads 4
```
**Step 7: Get allele-specific read counts with GATK** This step will filter a phased and genotyped vcf for heterozygous SNV to perform allele-specific counting in a coordinate-sorted and indexed bam file after WASP realignment. There are multiple options for count table outputs. The --pseudobulk option will group all barcodes together to perform allele-specific counting. This is analogous to bulk RNA-seq. The --celltype_counts option will use the barcode annotations to split the bam into cell-type-specific bam files before performing allele-specific counting. The --single_cell_counts option will split the input into individual single cell bam files and perform allele-specific counting. If your input vcf is phased you can select the --isphased option to add a phased genotype to the count tables.
```
Usage: step7_gatk_alleleCount.sh [-viognmlCcspt]
   -v  | --inputvcf           STR   path/to/input.vcf.gz eg. [vcfdir/funcotation/sample_1.pass.joint.hcphase.funco.vcf.gz]
   -i  | --inputbam           STR   path/to/wasp.bam eg. [project/wasp_rna/sample_1.phase.wasp.bam]
   -o  | --outputdir          STR   path/to/output directory eg. [project/wasp_rna/counts]
   -b  | --barcodes           STR   path/to/barcodes.csv eg. [barcodes/rna_barcodes.csv]
   -g  | --genotype           STR   genotype: [rna] [atac] [joint]
   -n  | --library_id         STR   library_id: eg. [sample_1]
   -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
   -l  | --interval           STR   optional: count a specified interval eg. [chr22]
   -C  | --pseudobulk_counts        allele-specific counts with all cells grouped together
   -c  | --celltype_counts          allele-specific counts after grouping cells by barcode annotation
   -s  | --single_cell_counts       single cell allele-specific counts for provided barcodes
   -p  | --isphased                 optional: input vcf is phased. Default=[false]
   -t  | --threads            INT   number of threads. Default=[1]
   -h  | --help                     show usage
```
**Launch 🌶️SALSA container**
```
SCRATCH1=/g/scratch
project=/g/salsa
reference=/g/reference
docker run \
--workdir $HOME \
-v $HOME:$HOME \
-v $reference/refdata-gex-GRCh38-2020-A:$HOME/rna_ref \
-v $project/vcf_output:$HOME/vcfdir \
-v $project/SALSA:$HOME/SALSA \
-v $project/barcodes:$HOME/barcodes \
-v $project:$HOME/project \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/g/scratch" \
--rm -it p4rkerw/salsa:latest
```
**Get phased allele-specific counts**
```
# runtime ~11min
bash SALSA/step7_gatk_alleleCount.sh \
--inputvcf vcfdir/funcotation/pbmc.pass.joint.chr22hcphase.funco.vcf.gz \
--inputbam project/wasp_rna/pbmc.hcphase.chr22wasp.bam \
--outputdir project/wasp_rna/counts \
--barcodes barcodes/rna_barcodes.csv \
--genotype rna_genotype \
--library_id pbmc_1k \
--modality rna \
--pseudobulk_counts \
--single_cell_counts \
--celltype_counts \
--interval chr22 \
--isphased \
--threads 10
```
**Inspect a single cell count table**
```
head -n5 project/wasp_rna/counts/pbmc_1k.chr22counts.single_cell.table
#contig  position        variantID       refAllele       altAllele       refCount        altCount        totalCount      lowMAPQDepth    lowBaseQDepth   rawDepth     otherBases # improperPairs    barcode
# chr22   50223608        .       G       A       1       0       1       0       0       1       0       0       AATCACGAGGAACTCG-1
# chr22   37282510        .       A       G       0       1       1       0       0       1       0       0       AATCACGCACTACCGG-1
# chr22   36281868        .       A       G       0       1       1       0       0       1       0       0       AATCACGGTATAGGAT-1
# chr22   36767036        .       C       T       0       1       1       0       0       1       0       0       AATCACGGTATAGGAT-1












