# **🌶️SALSA! Single Cell ATAC Tutorial**

**🌶️SALSA** is a tool for generating and analyzing phased single cell allele-specific read counts from 10X Genomics cellranger datasets. For additional information, please consult the [10X Genomics website](https://www.10xgenomics.com/). The workflow runs in a publicly-available docker container with all the necessary dependencies for code execution. For this tutorial, we will download a dataset from the 10X Genomics website. To complete all of the steps in the tutorial, you will need to download a cellranger reference, GATK bundle resources, and 1000G phased reference (see below for more information). All of the tutorial outputs are included in the this repository or on the 10X Genomics website in case you want to skip a step or compare results. If you are analyzing your own data, you will need to run cellranger and annotate your barcodes before starting the workflow.

**Step 0: Download cellranger reference** If you don't already have a GRCh38 cellranger reference download one from the [10X Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). 10X Genomics routinely updates their references with each new cellranger build, but new references are often backwards-compatible. The cellranger reference in this tutorial is compatible with the tutorial dataset (which is aligned to GRCh38-2020-A). Feel free to download a different reference and/or [dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets) from the 10X Genomics collection; just make sure it's aligned to GRCh38 so it matches the GATK bundle resources and ucsc contig style. Alignment information can be found in the summary html files. The tutorial assumes that you downloaded a reference to /mnt/g/reference so make sure you change the path if you used a different directory. 
```
reference=/mnt/g/reference
wget -P $reference https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

# check md5
md5sum $reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz #2f12f6016197869e21e5559827002716

# unpack
tar -xvzf $reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
```

**(Optional) Step 0: Download tutorial fastq files** We will download a a single cell ATAC dataset obtained from 1k PBMCs from a healthy donor. This dataset uses the NextGem v1.1 chemistry. If you want to skip ahead to the 🌶️SALSA workflow, you can download the coordinate-sorted bam and index from the 10X Genomics website. See Step 1 for details. 
```
# URL to the dataset: https://support.10xgenomics.com/single-cell-atac/datasets/2.0.0/atac_pbmc_1k_nextgem
# create your salsa tutorial directory and download the fastq
project=/mnt/g/salsa
wget -P $project/tar https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_fastqs.tar

# check md5
md5sum $project/tar/atac_pbmc_1k_nextgem_fastqs.tar #2a7af90faea21b25b6549ccb727106cc

# unpack
tar -C $project -xvf $project/tar/atac_pbmc_1k_nextgem_fastqs.tar
```
**(Optional) Step 0: Align the tutorial dataset to the cellranger reference** If you want to skip ahead to the 🌶️SALSA workflow, you can download the coordinate-sorted bam and index from the 10X Genomics website. See Step 1 for details. 
```
# align and count with cellranger-atac (version 2.0)
# the output will appear in the working directory which should be $project
reference=/mnt/g/reference
cellranger-atac count \
--id=pbmc_1k \
--reference=$SCRATCH1/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--fastqs=$SCRATCH1/salsa/atac_pbmc_1k_nextgem_fastqs \
--sample=atac_pbmc_1k_nextgem \
--localcores 10
```
**Step 0: Download GATK resource bundle** The following files are required for genotyping with GATK using GRCh38 and can be found in the [GATK google cloud bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false) . Download these files to a folder called gatk in your reference directory. For additional information on GATK germline and RNA-seq short variant discovery check out their [website](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
```
# download to /mnt/g/reference/gatk
# files needed for GATK germline short variant discovery
resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz
resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz
resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz
resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
resources_broad_hg38_v0_wgs_calling_regions.hg38
resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz
```

**Step 0: Pull 🌶️SALSA container** 
```
docker pull p4rkerw/salsa:latest
```
The 🌶️SALSA docker container is built on [broadinstitute/gatk:4.2.0.0](https://hub.docker.com/r/broadinstitute/gatk) with additional dependencies pre-installed (see below). Docker is a set of platform as a service (PaaS) products that use OS-level virtualization to deliver software in packages called containers, which can be used to enhance reproducibility. You can find the the 🌶️SALSA Dockerfile [here](https://github.com/p4rkerw/SALSA/blob/main/docker/Dockerfile).
```
# container dependencies
GATK 4.2.0.0
bwa 0.7.17
STAR 2.5.1b
bcftools 1.9 
pysam 0.15.3
shapeit 4.2
WASP 0.3.4
```

**Step 0: Clone 🌶️SALSA github repository** The repository is cloned to the $project directory, which is the same directory that the tutorial dataset was downloaded to. The tutorial assumes that the path to your project directory is /mnt/g/salsa so make sure to change the path if you chose a different directory.
```
project=/mnt/g/salsa
git -C $project clone https://github.com/p4rkerw/SALSA
```

**Launch 🌶️SALSA container** : Mount the required volumes in an interactive session. The $SCRATCH1 variable designates a temporary file directory. The tutorial assumes that your scratch directory is located at /mnt/g/scratch so make sure to change the path if you're using a different directory. The --cpus flag is included to limit the number of threads that each Java process can access within the container. Importantly, it does not limit the total number of threads that can be used by GATK. GATK runs on OpenJDK 1.8, which was the earliest version of Java to have container support. When a docker container is launched, each Java process can access all of the available threads within the container. Unfortunately, Java processes in docker are not aware of other Java processes running in parallel. As a result, GATK may exhaust your system resources when analyzing multiple genomic intervals. You probably wont need to set the --cpus flag in a in a high throughput computational environment (like a cluster), but its absence may cause Java to crash on a workstation. The --cpus flag is only needed for the genotyping step and should be left out in subsequent steps. You can read more about docker and Java [here](https://blog.softwaremill.com/docker-support-in-new-java-8-finally-fd595df0ca54)     
```
SCRATCH1=/mnt/g/scratch
project=/mnt/g/salsa
reference=/mnt/g/reference
docker run \
--cpus 2 \
--workdir $HOME \
-v $HOME:$HOME \
-v $project:$HOME/project \
-v $reference:$HOME/reference \
-v $project/SALSA:$HOME/SALSA \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/mnt/g/scratch" \
--rm -it p4rkerw/salsa:latest
```

**Step 0: Index and validate the GATK resource bundle files** We will use GATK to index the resource bundle files. In the process of indexing the vcfs, GATK will also validate them for proper formatting. 
```
# check the md5sum of the gatk resource files to ensure download integrity
md5sum reference/gatk/resources_broad_hg38_v0*
# b2979b47800b59b41920bf5432c4b2a0  reference/gatk/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf
# f7e1ef5c1830bfb33675b9c7cbaa4868  reference/gatk/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
# 14cc588a271951ac1806f9be895fb51f  reference/gatk/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf
# 2e02696032dcfe95ff0324f4a13508e3  reference/gatk/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf
# d05ac6b9a247a21ce0030c7494194da9  reference/gatk/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf
# 1790be6605825971526fff7cb3232764  reference/gatk/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list

# ensure the GATK resource bundle vcf are compressed prior to indexing. Indexing the files with GATK has the added benefit 
# stringent vcf format validation. Processing the dbsnp resource will take a few minutes...
for vcf in $(ls reference/gatk/*.vcf); do
  bcftools view -Oz $vcf > $vcf.gz
  gatk IndexFeatureFile -I $vcf.gz
  rm $vcf
done  

# check contents of reference/gatk folder
ls -1a reference/gatk
# resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz
# resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz
# resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
# resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz
# resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
# resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
# resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
# resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz
# resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz.tbi
# resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list

```

**Step 1: Genotype a single cell gene expression dataset** The tutorial workflow is based on the GATK germline short variant discovery pipeline. Additional info can be found on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) . To explore additional GATK options type 'gatk --list' into the terminal. If you want to skip ahead move the [genotyped vcf](https://github.com/p4rkerw/SALSA/blob/main/Tutorials/single_cell_atac/pbmc.atac.chr22.vcf.gz) and its index to the volume mounted to project/atac_genotype and proceed to the next step.
```
Usage: step1_gatk_genotype.sh [-inrgdomlVt]
  -i  | --inputbam           STR   path/to/input.bam eg. [rna_counts/sample_1/outs/possorted*.bam]
  -n  | --library_id         STR   library_id: eg. [sample_1]
  -r  | --reference          STR   path/to/cellranger_ref eg. [reference/refdata-gex-GRCh38-2020-A]
  -g  | --gatk_bundle        STR   path/to/gatk_bundle eg. [reference/gatk]
  -d  | --outputdir          STR   output directory name eg. [project/rna_genotype]
  -o  | --outputvcf          STR   name of output vcf eg. [sample_1.rna.vcf.gz]
  -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
  -l  | --interval           STR   optional: genotype a single chromosome eg. [chr22]
  -V  | --verbose                  optional: stream GATK output to terminal. Default=[false]
  -t  | --threads            INT   number of threads. Default=[1]
  -h  | --help                     show usage

```

**Genotype an ATAC sample with 🌶️SALSA** There are two options in the tutorial workflow depending on whether you did the cellranger alignment or skipped ahead. You only need to do one. 
```
# Option 1: if you followed the cellranger alignment step...
bash SALSA/step1_gatk_genotype.sh \
--inputbam project/pbmc_1k/outs/possorted_bam.bam \
--library_id pbmc_1k \
--reference reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--gatk_bundle reference/gatk \
--outputdir project/atac_genotype \
--outputvcf pbmc.atac.chr22.vcf.gz \
--interval chr22 \
--modality atac \
--threads 4

# Option 2: if you skipped the cellranger alignment step...
# download coordinate-sorted bam and index
wget -P project/pbmc_1k https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_possorted_bam.bam
wget -P project/pbmc_1k https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_possorted_bam.bam.bai

# genotype the downloaded bam
bash SALSA/step1_gatk_genotype.sh \
--inputbam project/pbmc_1k/atac_pbmc_1k_nextgem_possorted_bam.bam \
--library_id pbmc_1k \
--reference reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--gatk_bundle reference/gatk \
--outputdir project/atac_genotype \
--outputvcf pbmc.atac.chr22.vcf.gz \
--interval chr22 \
--modality atac \
--threads 4
```
**Inspect the genotyped vcf** Note how the GT field for the majority of variants has a "/" character, which indicates that the variants are not phased. The variants with a "|" are physically-phased by GATK HaplotypeCaller.
```
bcftools query -f '[%CHROM,%POS,%REF,%ALT,%GT,%FILTER\n]' project/atac_genotype/pbmc.atac.chr22.vcf.gz | head -n5
# chr22,10673696,A,G,1/1,PASS
# chr22,10713908,G,A,1|1,PASS
# chr22,10713916,G,A,1|1,PASS
# chr22,10714219,G,A,0/1,PASS
# chr22,10714597,G,C,1/1,PASS
```
**(Not required for tutorial) Step 2: Merge genotypes from the same patient** If you genotyped a paired single cell gene expression and ATAC dataset from a split sample (or a single cell Multiome) you can merge genotypes into a single vcf. If you're following the tutorial, you can skip this step.
```
Usage: step2_merge_geno.sh [-nabdit]
  -n  | --library_id         STR   library_id: eg. [sample_1]
  -a  | --vcfone             STR   path/to/first.vcf.gz eg. [project/atac_genotype/sample_1.atac.vcf.gz]. Prioritize annotations from first vcf if there is overlap.
  -b  | --vcftwo             STR   path/to/second.vcf.gz eg. [project/rna_genotype/sample_1.rna.vcf.gz]
  -d  | --outputdir          STR   output directory [project/joint_genotype]
  -o  | --outputvcf          STR   name of output vcf eg. [sample_1.pass.joint.vcf.gz]
  -i  | --includefiltered          optional: include filtered variants in output vcf. Default=[false]
  -t  | --threads            INT   number of threads. Default=[1]
  -h  | --help                     show usage
```

**Launch 🌶️SALSA container** : Note how the --cpus flag is not needed for the rest of the workflow    
```
SCRATCH1=/mnt/g/scratch
project=/mnt/g/salsa
reference=/mnt/g/reference
docker run \
--workdir $HOME \
-v $HOME:$HOME \
-v $project:$HOME/project \
-v $reference:$HOME/reference \
-v $project/SALSA:$HOME/SALSA \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/mnt/g/scratch" \
--rm -it p4rkerw/salsa:latest
```

**(Recommended) Step 3: Phase genotype with 🌶️SALSA** If you want to perform your analysis with phased genotypes you will need a phased reference. We will use [shapeit4](https://github.com/odelaneau/shapeit4) to phase our variants. To explore additional phasing options type 'shapeit4.2' into your terminal. Variant phasing increases the performance of the WASP variant-realignment and downstream analysis steps. Download the 1000G phased reference files for SNV only or SNV_and_INDEL from ftp.1000genomes.ebi.ac.uk . If you are analyzing the tutorial ATAC dataset select the SNV reference. If you want to skip ahead download the [phased vcf](https://github.com/p4rkerw/SALSA/blob/main/Tutorials/single_cell_atac/pbmc.pass.atac.chr22hcphase.vcf.gz) and its index to the volume mounted to project/phasing and proceed to the next step.

a) SNV only: [/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/) </br>
b) SNV_and_INDEL: [/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/)
```
Usage: step3_phase_vcf.sh [-nvdorlpsitrh]
  -n  | --library_id         STR   library_id: eg. [sample_1]
  -v  | --inputvcf           STR   path/to/input.vcf.gz eg. [project/joint_genotype/sample_1.pass.joint.vcf.gz]
  -d  | --outputdir          STR   output directory name eg. [project/phasing]
  -o  | --outputvcf          STR   name of output vcf eg. [sample_1.pass.joint.phase.vcf.gz]
  -r  | --phasingref         STR   path/to/1000G reference eg. [reference/phasing/biallelic_SNV]
  -l  | --interval           STR   optional: phase a single chromosome eg. [chr22]
  -p  | --hcphase                  optional: recover haplotypecaller physical phasing variants that are not in shapeit reference. Default=[false]
  -s  | --snvonly                  use the biallelic_SNV reference for phasing
  -i  | --snvindel                 use the biallelic_SNV_and_INDEL reference for phasing
  -r  | --reproduce                optional: run shapeit with a single thread for reproducibility. Default=[false]
  -t  | --threads            INT   number of threads. Default=[1]
  -h  | --help                     show usage
```

You will eventually need to download the vcf for every chromosome, but for the purposes of the tutorial just download the SNV reference for chr22. The reference may take awhile to download so feel free to move ahead to the next step using the [phased vcf](https://github.com/p4rkerw/SALSA/blob/main/Tutorials/single_cell_atac/pbmc.pass.atac.chr22hcphase.vcf.gz) in the repository. You can come back to step 3 when the download finishes.
```
# download chr22 SNV phasing reference
url=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget -P reference/phasing/biallelic_SNV $url
```

**Rename the reference contigs** The 1000G vcf reference files do not have the same contig style as the cellranger reference. You will need to update the 1000G contig style using bcftools. For the tutorial, we will only do chromosome 22. 
```
# runtime ~10min
# rename the contigs in the 1000G reference from 22 to chr22
for i in {1..22} X;do echo "${i} chr${i}";done > /tmp/rename_chrm.txt
inputvcf=ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
bcftools annotate reference/phasing/biallelic_SNV/$inputvcf --threads 4 --rename-chrs /tmp/rename_chrm.txt -Oz -o $SCRATCH1/$inputvcf
mv $SCRATCH1/$inputvcf reference/phasing/biallelic_SNV
bcftools index --threads 4 reference/phasing/biallelic_SNV/$inputvcf
```
**Phase an interval with 🌶️SALSA** For the tutorial, we will use multithreading to phase the variants, which means we will not set the --reproduce flag. As a result, there may be small differences between the phased vcf in the repository and your vcf. 
```
# runtime ~9min
bash SALSA/step3_phase_vcf.sh \
--library_id pbmc_1k \
--inputvcf project/atac_genotype/pbmc.atac.chr22.vcf.gz \
--outputdir project/phasing \
--outputvcf pbmc.pass.atac.chr22hcphase.vcf.gz \
--phasingref reference/phasing/biallelic_SNV \
--interval chr22 \
--hcphase \
--snvonly \
--threads 10
```
**Inspect the phased vcf** The GT field has a "|" character, indicating that the variants are now phased. We did not set the --reproduce flag so there may be small differences between these results and your results. 
```
bcftools query -f '[%CHROM,%POS,%REF,%ALT,%GT\n]' project/phasing/pbmc.pass.atac.chr22hcphase.vcf.gz | head -n5
# chr22,11030215,G,A,0|1
# chr22,11312610,C,T,1|1
# chr22,11355712,C,T,1|1
# chr22,11372501,A,G,1|1
# chr22,12602521,C,T,1|1
```
**(Optional) Step 4: Annotate vcf with GATK Funcotator** If you want to annotate your vcf with GENCODE and gnomAD you will need to download the [GATK Funcotator resource](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial). GATK routinely updates its resources so you may need to change the name of the folder in the tutorial to match the one you downloaded. gnomAD resources need to be enabled after download (see GATK instructions on their website). When the resources have been downloaded, move the dataSources folder into to the reference directory (eg. [reference/funcotator_dataSources.v1.6.20190124]). If you want to skip ahead while these files are downloading move the [funcotated vcf](https://github.com/p4rkerw/SALSA/blob/main/Tutorials/single_cell_atac/pbmc.pass.atac.chr22hcphase.funco.vcf.gz) and its index to the volume mounted to project/funcotation and proceed to the next step.
```
Usage: step4_gatk_anno_vcf.sh [-nvdoramfVth]
  -n  | --library_id         STR   library_id: eg. [sample_1]
  -v  | --inputvcf           STR   path/to/input.vcf.gz eg. [project/phasing/sample_1.pass.joint.hcphase.vcf.gz]
  -d  | --outputdir          STR   output directory name eg. [project/funcotation]
  -o  | --outputvcf          STR   name of output vcf eg. [sample_1.pass.joint.hcphase.funco.vcf.gz]
  -r  | --reference          STR   path/to/cellranger_ref eg. [reference/refdata-gex-GRCh38-2020-A]
  -a  | --output_table       STR   name of output funcotation csv eg. [sample_1.pass.joint.hcphase.formatted.csv]
  -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
  -f  | --funcotation        STR   path/to/funcotation directory eg. [reference/funcotator_dataSources.v1.6.20190124g]
  -V  | --verbose                  optional: stream GATK output to terminal. Default=[false]
  -t  | --threads            INT   number of threads. Default=[1]
  -h  | --help                     show usage
```

**Annotate a vcf with 🌶️SALSA**
```
# runtime ~1min
bash SALSA/step4_gatk_anno_vcf.sh \
--library_id pbmc_1k \
--inputvcf project/phasing/pbmc.pass.atac.chr22hcphase.vcf.gz \
--outputdir project/funcotation \
--outputvcf pbmc.pass.atac.chr22hcphase.funco.vcf.gz \
--reference reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--funcotation reference/funcotator_dataSources.v1.6.20190124g \
--output_table pbmc.pass.atac.chr22hcphase.formatted.csv \
--modality atac \
--threads 10
```
**Inspect the annotation table** GATK Funcotator provides a lot of annotation fields in the vcf INFO field and only a subset are included in this table. These two variants are located in an intergenic region.
```
head -n3 project/funcotation/pbmc.pass.atac.chr22hcphase.formatted.csv
# variant_id,CHROM,POS,REF,ALT,GT,FILTER,Gencode_27_variantClassification,Gencode_27_codonChange,gnomAD_exome_AF,gnomAD_genome_AF,Gencode_27_hugoSymbol
# chr22_11030215_G_A,chr22,11030215,G,A,0|1,.,IGR,,,,Unknown
# chr22_11312610_C_T,chr22,11312610,C,T,1|1,.,IGR,,,,Unknown
```

**(Recommended) Step 5: Filter the bam file with cell barcodes** Use barcode celltype annotations to filter the coordinate-sorted cellranger bam using the CB tag. This step will speed up downstream analysis by eliminating barcodes that do not meet quality control. The barcode annotation file should have three columns where the first column is the barcode, the second column is the library_id, and the third column is the celltype annotation. For the purposes of the tutorial, we will only filter chr22.
```
Usage: step5_filterbam.sh [-nidolmbeth]
-n  | --library_id         STR   library_id: eg. [sample_1]
-i  | --inputbam           STR   path/to/input.bam eg. [project/sample_1/outs/possorted_genome_bam.bam]
-d  | --outputdir          STR   output directory eg. [project/wasp_rna]
-o  | --outputbam          STR   filtered output bam eg. [sample_1.bcfilter.bam]
-l  | --interval           STR   optional: filter a single chromosome eg. [chr22]
-m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
-b  | --barcodes           STR   path/to/barcodes.csv with headers and two columns. First column is named "barcodes"
                                 second column is group "orig.ident" and third column is cell type eg. [project/barcodes/rna_barcodes.csv]
-e  | --validate                 validate the barcode-filtered bam file. Default=[false]
-t  | --threads            INT   number of threads. Default=[1]
-h  | --help                     show usage
```
**Download clustering analysis for tutorial dataset and create a barcode csv**
```
# download the barcode cluster annotation file from 10X Genomics
wget -P project/pbmc_1k https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_analysis.tar.gz
tar -C project/pbmc_1k -xvzf project/pbmc_1k/atac_pbmc_1k_nextgem_analysis.tar.gz

# use the cluster number as celltype and assign pbmc_1k as the library_id in final csv
cluster=project/pbmc_1k/analysis/clustering/graphclust/clusters.csv
mkdir project/barcodes
(echo "barcode,orig.ident,celltype"; awk 'BEGIN{FS=OFS=","} {if (NR!=1) print $1,"pbmc_1k",$2}' $cluster) |sed 's/-1//g'> project/barcodes/atac_barcodes.csv
```
**Inspect the barcode csv**
```
head -n3 project/barcodes/atac_barcodes.csv
# barcode,orig.ident,celltype
# AAACGAATCGCATAAC,pbmc_1k,6
# AAACGAATCTGTGTGA,pbmc_1k,1

```
**Filter cellranger bam with barcode csv**
```
# runtime ~1min
bash SALSA/step5_filterbam.sh \
--library_id pbmc_1k \
--inputbam project/pbmc_1k/outs/possorted_bam.bam \
--modality atac \
--interval chr22 \
--barcodes project/barcodes/atac_barcodes.csv \
--outputdir project/wasp_atac \
--outputbam pbmc.bcfilter.chr22.bam \
--threads 10
```
**Step 6: Perform variant-aware realignment with WASP** This step takes a genotyped vcf and performs variant-aware realignment on a coordinate-sorted and indexed bam file with WASP. WASP is a tool to perform unbiased allele-specific read mapping and you can read more about it [here](https://github.com/bmvdgeijn/WASP). For the purposes of the tutorial, we will only analyze chromosome 22. For ATAC analysis, we will build a bwa index of the cellranger reference at run time. 
```
Usage: step6_wasp.sh [-vbdogianlmpt]
  -v  | --inputvcf          STR   project/funcotation/sample_1.pass.joint.hcphase.funco.vcf.gz
  -b  | --inputbam          STR   path/to/input.bam eg. [project/wasp_rna/sample_1.bcfilter.bam]
  -d  | --outputdir         STR   name of output directory eg. [project/wasp_rna]
  -o  | --outputbam         STR   name of output wasp bam eg. [sample_1.phase.wasp.bam]
  -g  | --genotype          STR   genotype: [rna] [atac] [joint]
  -i  | --stargenome        STR   path/to/star genome index for STAR alignment eg. [reference/refdata-gex-GRCh38-2020-A/star]
  -a  | --atacref           STR   path/to/atac_reference for bwa alignment eg. [reference/refdata-cellranger-atac-GRCh38-1.2.0]
  -n  | --library_id        STR   library_id: eg. [sample_1]
  -m  | --modality          STR   modality: [rna] [atac]
  -l  | --interval          STR   optional: analyze a single chromosome eg. [chr22]
  -p  | --isphased                input vcf is phased. Default=[false]
  -t  | --threads           INT   number of threads. Default=[1]
  -h  | --help                    show usage
```

**Download the hg38 chromInfo file from UCSC**
```
reference=/mnt/g/reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz -O $reference/hg38_chromInfo.txt.gz
```

**Run WASP on the barcode-filtered bam**
```
# runtime ~2min
bash SALSA/step6_wasp.sh \
--inputvcf project/funcotation/pbmc.pass.atac.chr22hcphase.funco.vcf.gz \
--inputbam project/wasp_atac/pbmc.bcfilter.chr22.bam \
--outputdir project/wasp_atac \
--outputbam pbmc.hcphase.chr22wasp.bam \
--atacref reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--genotype atac \
--library_id pbmc_1k \
--modality atac \
--isphased \
--interval chr22 \
--threads 10
```
**Step 7: Get allele-specific read counts with 🌶️SALSA** This step will filter a phased and genotyped vcf for heterozygous SNV to perform allele-specific counting in a coordinate-sorted and indexed bam file after WASP realignment. There are multiple options for count table outputs. The --pseudobulk option will group all barcodes together to perform allele-specific counting. This is analogous to bulk ATAC-seq. The --celltype_counts option will use the barcode annotations to split the bam into cell-type-specific bam files before performing allele-specific counting. The --single_cell_counts option will split the input into individual single cell bam files and perform allele-specific counting. If your input vcf is phased you can select the --isphased option to add a phased genotype to the count tables.
```
Usage: step7_gatk_alleleCount.sh [-viognmrlCcspt]
  -v  | --inputvcf           STR   path/to/input.vcf.gz eg. [project/funcotation/sample_1.pass.joint.hcphase.funco.vcf.gz]
  -i  | --inputbam           STR   path/to/wasp.bam eg. [project/wasp_rna/sample_1.phase.wasp.bam]
  -b  | --barcodes           STR   path/to/barcodes.csv eg. [project/barcodes/rna_barcodes.csv]
  -g  | --genotype           STR   genotype: [rna] [atac] [joint]
  -n  | --library_id         STR   library_id: eg. [sample_1]
  -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
  -r  | --reference          STR   path/to/cellranger_ref eg. [reference/refdata-gex-GRCh38-2020-A]
  -l  | --interval           STR   optional: count a specified chromosome eg. [chr22]
  -C  | --pseudobulk_counts        allele-specific counts with all cells grouped together
  -c  | --celltype_counts          allele-specific counts after grouping cells by barcode celltype annotation
  -s  | --single_cell_counts       single cell allele-specific counts for provided barcodes     
  -p  | --isphased                 optional: input vcf is phased. Default=[false]
  -V  | --verbose                  optional: stream GATK output to terminal. Default=[false]
  -t  | --threads            INT   number of threads. Default=[1]
  -h  | --help                     show usage
```
**Get phased allele-specific counts with 🌶️SALSA**
```
# runtime ~6min
bash SALSA/step7_gatk_alleleCount.sh \
--inputvcf project/funcotation/pbmc.pass.atac.chr22hcphase.funco.vcf.gz \
--inputbam project/wasp_atac/pbmc.hcphase.chr22wasp.bam \
--outputdir project/wasp_atac/counts \
--barcodes project/barcodes/atac_barcodes.csv \
--genotype atac_genotype \
--library_id pbmc_1k \
--modality atac \
--reference reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--pseudobulk_counts \
--single_cell_counts \
--celltype_counts \
--interval chr22 \
--isphased \
--threads 10
```
**Inspect a single cell count table**
```
head -n5 project/wasp_atac/counts/pbmc_1k.chr22counts.single_cell.table
# contig  position        variantID       refAllele       altAllele       refCount        altCount        totalCount      lowMAPQDepth    lowBaseQDepth   rawDepth        otherBases      improperPairs       barcode
# chr22   23751238        .       A       T       0       1       1       0       0       1       0       0       AAACGAATCGCATAAC-1
# chr22   23914067        .       C       T       1       0       1       0       0       1       0       0       AAACGAATCGCATAAC-1
# chr22   36507049        .       C       T       0       1       1       0       0       1       0       0       AAACGAATCGCATAAC-1
# chr22   43187819        .       A       C       0       1       1       0       0       1       0       0       AAACGAATCGCATAAC-1

