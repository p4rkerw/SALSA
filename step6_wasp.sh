#!/bin/bash
# this script will take a vcf and a barcode-filtered cellranger bam file and perform variant aware realignment with WASP

# Set some default values:
inputvcf=""
genotype=""
stargdir=""
atacref=""
library_id=""
modality=""
isphased=false
threads=1
r2only=false
interval=""

function usage {
cat << "EOF"


        /|      (                (      (              
     .-((--.     )\ )     (       )\ )   )\ )     (    
    ( '`^'; )   (()/(     )\     (()/(  (()/(     )\   
    `;#    |     /(_)) ((((_)(    /(_))  /(_)) ((((_)( 
     \#    |    (_))    )\ _ )\  (_))   (_))    )\ _ )\ 
      \#   \    / __|   (_)_\(_) | |    / __|   (_)_\(_) 
       '-.  )   \__ \    / _ \   | |__  \__ \    / _ \   
          \(    |___/   /_/ \_\  |____| |___/   /_/ \_\ 
           `

Single Cell Allele Specific Analysis
Author: Parker C. Wilson MD, PhD
Contact: parkerw@wustl.edu
Version: 1.0

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

EOF
exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

PARSED_ARGUMENTS=$(getopt -a -n step6_wasp.sh \
-o v:b:d:o:g:i:a:n:m:l:pt:h \
--long inputvcf:,inputbam:,outputdir:,outputbam:,genotype:,stargenome:,atacref:,library_id:,modality:,interval:,isphased,threads:,help -- "$@")

echo "PARSED_ARGUMENTS are $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -v | --inputvcf)          inputvcf=$2      ; shift 2 ;;
    -b | --inputbam)          inputbam=$2      ; shift 2 ;;
    -d | --outputdir)         outputdir=$2     ; shift 2 ;;
    -o | --outputbam)         outputbam=$2     ; shift 2 ;;
    -g | --genotype)          genotype=$2      ; shift 2 ;;
    -i | --stargenome)        stargdir=$2      ; shift 2 ;;
    -a | --atacref)           atacref=$2       ; shift 2 ;;
    -n | --library_id)        library_id=$2    ; shift 2 ;;
    -m | --modality)          modality=$2      ; shift 2 ;;
    -l | --interval)          interval=$2      ; shift 2 ;;
    -p | --isphased)          isphased=true    ; shift 1 ;;
    -t | --threads)           threads=$2       ; shift 2 ;;
    -h | --help)              usage ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

echo "library_id   : $library_id"
echo "inputvcf     : $inputvcf"
echo "inputbam     : $inputbam"
echo "outputdir    : $outputdir"
echo "outputbam    : $outputbam"
echo "stargenome   : $stargdir"
echo "atacref      : $atacref"
echo "modality     : $modality"
echo "interval     : $interval"
echo "isphased     : $isphased"
echo "threads      : $threads"
echo "Parameters remaining are: $@"

# checking input files
if [ ! -f $inputbam ]; then { echo "Input bam file not found"; exit 1; }; fi
if [ ! -f $inputvcf ]; then { echo "Input vcf file not found"; exit 1; }; fi

# TODO: check for reference

# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

contigdir=/tmp/vcf/$library_id
wasp=/opt/WASP
workdir=$SCRATCH1/wasp_${modality}/${genotype}_genotype/find_intersecting_snps_${library_id}
mkdir -p $workdir 
mkdir -p $outputdir 2> /dev/null
mkdir -p $contigdir

# following command returns 1000 for paired-end library and 0 for single-end
read_config=$((samtools view -H $inputbam ; samtools view $inputbam |head -n1000) | samtools view -c -f 1)
echo "Sampling first 1000 reads of input bam to determine read confguration"
if [ $read_config -eq 0 ]; then
  echo "Detected single-end bam"
  r2only=true
elif [ $read_config -eq 1000 ]; then
  echo "Detected paired-end bam"
  is_paired_end="--is_paired_end"
fi

# index the input bam if no index is detected
if [ ! -f $inputbam.bai ]; then
  samtools index -@ $threads $inputbam
fi

# create a basename for temporary files. find_intersecting_snps.py does not allow naming of outputs and outputs files with a basename
bn=$(basename $inputbam .bam)

# filter by selected interval if specified
bcftools index --threads $threads --tbi $inputvcf 2> /dev/null
if [ $interval ]; then
  echo "Selected interval is ${interval}"
  workdir=$workdir/$interval
  mkdir -p $workdir
  # filter input bam for selected interval
  echo "Filtering bam for selected interval $interval"
  samtools view -@ $threads -bS $inputbam $interval|samtools sort -@ $threads -T $workdir > $workdir/$bn.bam
  inputbam=$workdir/$bn.bam
  bn=$(basename $inputbam .bam)
  samtools index -@ $threads $inputbam
  contigs=($interval)
else
  # generate array of unique contigs
  echo "Creating contig files from vcf"
  contigs=($(bcftools query -f'[%CHROM\n]' $inputvcf|grep chr|sort|uniq))
fi


if [ $isphased = "true" ]; then
  # set output file name to indicate that it's phased
  echo "Input vcf is phased"
  # chrominfo from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz
  # first separate input vcf into contigs as required by snp2h5 utility
  rm /tmp/haplotypes.h5 /tmp/snp_index.h5 /tmp/snp_tab.h5 /tmp/haplotype.chr*.vcf.gz 2> /dev/null
  for contig in ${contigs[@]}; do
    echo "Subsetting vcf for $contig"
    bcftools view -Oz $inputvcf $contig -o /tmp/haplotype.$contig.vcf.gz
  done
  # create snp hdf5 file
  echo "Generating snp hdf5"
  $wasp/snp2h5/snp2h5 \
    --chrom reference/hg38_chromInfo.txt.gz \
    --format vcf \
    --haplotype /tmp/haplotypes.h5 \
    --snp_index /tmp/snp_index.h5 \
    --snp_tab   /tmp/snp_tab.h5 \
    /tmp/haplotype.chr*.vcf.gz
  # the is_paired_end variable expands to --is_paired_end for paired-end reads
  echo "Running WASP and writing to $workdir"
  python $wasp/mapping/find_intersecting_snps.py \
    ${is_paired_end} \
    --is_sorted \
    --output_dir $workdir \
    --snp_index /tmp/snp_index.h5 \
    --snp_tab /tmp/snp_tab.h5 \
    --haplotype /tmp/haplotypes.h5 \
    $inputbam
fi


# for unphased vcf use a text-based snv file as recommended by WASP
if [ $isphased = "false" ]; then
  # set output file name to indicate that it was not phased
  # split the vcf into separate files by contig and print ref and alt for each variant
  contigdir=/tmp/vcf 
  mkdir $contigdir 2> /dev/null
  for contig in $contigs; do
    echo "Generating snv contig file for $contig"
    output_file=$contigdir/$contig.snps.txt.gz
    # get SNPs from VCF files:
    bcftools view -H $inputvcf|awk -v a=$contig '{if($1 == a) {print $2,$4,$5}}'|gzip > $output_file
  done
  # the is_paired_end variable expands to --is_paired_end for paired-end reads
  echo "Running WASP and writing to $workdir"
  python $wasp/mapping/find_intersecting_snps.py \
    ${is_paired_end} \
    --is_sorted \
    --output_dir $workdir \
    --snp_dir $contigdir \
    $inputbam
fi


################SCRNA SPECIFIC WORKFLOW################
# create a STAR index if not specified
if [ $modality = "rna" ]; then
  if [ ! -f $stargdir/genomeParameters.txt ]; then
    echo "Generating STAR index and putting in $stargdir directory"
    STAR \
      --runMode genomeGenerate \
      --runThreadN $threads \
      --genomeDir $stargdir \
      --genomeFastaFiles rna_ref/fasta/genome.fa \
      --sjdbGTFfile rna_ref/genes/genes.gtf
  fi
  ### use STAR to remap single end reads 
  if [ $r2only = "true" ]; then
    echo "Running STAR with single end reads"
    STAR \
      --genomeDir $stargdir \
      --runThreadN $threads \
      --readFilesIn <(gunzip -c $workdir/$bn.remap.fq.gz) \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix $workdir/$bn.
    mv $workdir/$bn.Aligned.sortedByCoord.out.bam $workdir/$bn.sorted.realigned.bam
    samtools index -@ $threads $workdir/$bn.sorted.realigned.bam
  elif [ $r2only = "false" ]; then
    echo "Running STAR with paired end reads"
    ### use STAR to remap paired end reads 
    STAR \
      --genomeDir $stargdir \
      --runThreadN $threads \
      --readFilesIn <(gunzip -c $workdir/$bn.remap.fq1.gz) <(gunzip -c $workdir/$bn.remap.fq2.gz) \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix $workdir/$bn.
    mv $workdir/$bn.Aligned.sortedByCoord.out.bam $workdir/$bn.sorted.realigned.bam
    samtools index -@ $threads $workdir/$bn.sorted.realigned.bam
  fi
fi


###########SNATAC SPECIFIC WORKFLOW###########
if [ $modality = "atac" ]; then
  # check for bwa mem index
  if [ ! -f $atacref/fasta/genome.dict ]; then
    echo "Generating BWA index and putting in atac_ref directory"
    bwa index $atacref/fasta/genome.fa 
  fi  
  # use bwa to remap reads
  echo "Realigning reads with BWA" 
  bwa mem -t $threads $atacref/fasta/genome.fa $workdir/$bn.remap.fq1.gz $workdir/$bn.remap.fq2.gz > $workdir/$bn.realigned.sam
  # sort the realigned bam file
  samtools view -bS $workdir/$bn.realigned.sam > $workdir/$bn.realigned.bam
  samtools sort -@ $threads -o $workdir/$bn.sorted.realigned.bam $workdir/$bn.realigned.bam
  samtools index -@ $threads $workdir/$bn.sorted.realigned.bam
fi


#################FILTER REMAPPED READS#################
# filter the remapped/realigned reads and output as keep.bam
python $wasp/mapping/filter_remapped_reads.py \
  $workdir/$bn.to.remap.bam \
  $workdir/$bn.sorted.realigned.bam \
  $workdir/keep.bam

###################MERGE AND SORT###################
echo "Merging realigned bam file"
# merge the output of find_intersecting_snps $bn.keep.bam and filter_remapped_reads keep.bam
samtools merge -f --threads $threads $workdir/$bn.keep.merge.bam \
  $workdir/$bn.keep.bam \
  $workdir/keep.bam
echo "Sorting output: $outputbam"  
samtools sort -@ $threads -T $workdir -o $outputdir/$outputbam \
  $workdir/$bn.keep.merge.bam
echo "Indexing output: $outputbam"
samtools index $outputdir/$outputbam
echo "Writing to directory: $outputdir"          

# NOTE: the WASP tool does not account for UMI / barcodes when it removes duplicates
# SOLUTION: retain cellranger duplicate markings and ASEReadCounter will automatically filter them OR
# split the bam by barcode and filter individual bams

# cleanup
rm -rf $workdir

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

echo -e "\033[35;40mWASP completed in $(format_time $SECONDS)\033[0m"
