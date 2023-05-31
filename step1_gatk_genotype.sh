#!/bin/bash
# this script will generate a genotyped vcf using the input bam file
# for info to implement a scatter workflow: https://www.ibm.com/downloads/cas/ZJQD0QAL

# exit on error
set -e

# Set some default values:
interval=""
threads=1
outputbam=""
verbose=false
exit_status=0

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

Usage: step1_gatk_genotype.sh [-inrgdomlbVt]
  -i  | --inputbam           STR   path/to/input.bam eg. [project/sample_1/outs/possorted*.bam]
  -n  | --library_id         STR   library_id: eg. [sample_1]
  -r  | --reference          STR   path/to/cellranger_ref eg. [reference/refdata-gex-GRCh38-2020-A]
  -g  | --gatk_bundle        STR   path/to/gatk_bundle eg. [reference/gatk]
  -d  | --outputdir          STR   output directory name eg. [project/rna_genotype]
  -o  | --outputvcf          STR   name of output vcf eg. [sample_1.rna.vcf.gz]
  -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
  -l  | --interval           STR   optional: genotype a single chromosome eg. [chr22]
  -b  | --outputbam          STR   optional: save bqsr bam to output dir. eg. [sample_1.rna.chr22.bam]
  -V  | --verbose                  optional: stream GATK output to terminal. Default=[false]
  -t  | --threads            INT   number of threads. Default=[1]
  -h  | --help                     show usage

EOF
exit 1
}

if [[ ${#} -eq 0 ]]; then usage; fi

PARSED_ARGUMENTS=$(getopt -a -n step1_gatk_genotype.sh \
-o i:n:r:g:d:o:m:l:b:Vt:h \
--long inputbam:,library_id:,reference:,gatk_bundle:,outputdir:,outputvcf:,modality:,interval:,outputbam:,verbose,threads:,help -- "$@")

echo "PARSED_ARGUMENTS are $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -i | --inputbam)            inputbam=$2                     ; shift 2 ;;
    -n | --library_id)          library_id=$2                   ; shift 2 ;;
    -r | --reference)           reference=$2                    ; shift 2 ;;
    -g | --gatk_bundle)         gatk_bundle=$2                  ; shift 2 ;;
    -d | --outputdir)           outputdir=$2                    ; shift 2 ;;
    -o | --outputvcf)           outputvcf=$2                    ; shift 2 ;;
    -m | --modality)            modality=$2                     ; shift 2 ;;
    -l | --interval)            interval=$2                     ; shift 2 ;;
    -b | --outputbam)           outputbam=$2                    ; shift 2 ;;
    -V | --verbose)             verbose=true                    ; shift 1 ;;
    -t | --threads)             threads=$2                      ; shift 2 ;;
    -h | --help)                usage ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

echo "inputbam                  : $inputbam"
echo "library_id                : $library_id"
echo "reference                 : $reference"
echo "gatk_bundle               : $gatk_bundle"
echo "outputdir                 : $outputdir"
echo "outputvcf                 : $outputvcf"
echo "modality                  : $modality"
echo "interval                  : $interval"
echo "outputbam                 : $outputbam"
echo "verbose                   : $verbose"
echo "threads                   : $threads"
echo "Parameters remaining are  : $@"

#####################################################
#####################################################
function gatk_germline_short_variant_scatter_gather {
# input bam is first positional argument $1
# bundle files can be obtained from gatk resource bundle on google cloud
# split selected intervals across the number of threads
# generate base recalibration table
  echo "Running BaseRecalibrator on scattered intervals for $interval"
  for scatter_interval in ${scatter_intervals[@]}; do
    gatk --java-options "-Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=4" BaseRecalibrator \
      -I $1 \
      -R $reference/fasta/genome.fa \
      --known-sites $gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
      --known-sites $gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
      --known-sites $gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
      --known-sites $gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
      -L /tmp/interval_files_folder/$scatter_interval \
      -O $intervaldir/recal_data.$scatter_interval.table \
      --tmp-dir $tmpdir \
      --verbosity INFO >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33mBaseRecalibrator failed on $interval. Check $outputlog for additional info\033[0m"; echo 1 > /tmp/exit_status.txt; } &
  done
  # exit with 1 if interval failed
  wait
  exit_status=$(head -n1 /tmp/exit_status.txt)
  if [ $exit_status -eq 1 ]; then return 1; fi

  # apply base quality score recalibration
  echo "Running ApplyBQSR on scattered intervals for $interval"
  for scatter_interval in ${scatter_intervals[@]}; do
    gatk --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap" ApplyBQSR \
      -I $1 \
      -R $reference/fasta/genome.fa \
      --bqsr-recal-file $intervaldir/recal_data.$scatter_interval.table \
      -L /tmp/interval_files_folder/$scatter_interval \
      -O $intervaldir/bqsr.$scatter_interval.bam \
      --tmp-dir $tmpdir \
      --verbosity INFO >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33mApplyBQSR failed on $interval. Check $outputlog for additional info\033[0m"; echo 1 > /tmp/exit_status.txt; } &
  done
  # exit with 1 if interval failed
  wait
  exit_status=$(head -n1 /tmp/exit_status.txt)
  if [ $exit_status -eq 1 ]; then return 1; fi

  # make sure bam outputs are indexed
  # and create list of bqsr interval bam
  if [ -f /tmp/ibam.list ]; then rm /tmp/ibam.list; fi
  for scatter_interval in ${scatter_intervals[@]}; do
    echo "$intervaldir/bqsr.$scatter_interval.bam" >> /tmp/ibam.list
    if [ ! -f $intervaldir/bqsr.$scatter_interval.bam.bai ]; then
      samtools index -@ $threads $intervaldir/bqsr.$scatter_interval.bam \
        || { echo -e "\033[0;33msamtools index failed on $scatter_interval\033[0m"; echo 1 > /tmp/exit_status.txt; }
    fi
  done
  exit_status=$(head -n1 /tmp/exit_status.txt)
  if [ $exit_status -eq 1 ]; then return 1; fi

  # gather the analysis-ready bam intervals and index
  if [ $outputbam ]; then
    echo "Gathering scattered bqsr bams for $interval"
    gatk GatherBamFiles \
      -I /tmp/ibam.list \
      -O $workdir/bqsr.$interval.bam >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33mGatherBamFiles failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }
    echo "Sorting gathered bam file"
    samtools sort -@ $threads $workdir/bqsr.$interval.bam -o $workdir/sorted.bqsr.$interval.bam >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33msamtools sort failed on $interval. Check $outputlog file for additional info\033[0m"; return 1; }
    samtools index -@ $threads $workdir/sorted.bqsr.$interval.bam >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33msamtools index failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }
    rm $workdir/bqsr.$interval.bam
  fi

  # remove vcf list if session has been used multiple times
  if [ -f /tmp/gvcf.list ]; then rm /tmp/gvcf.list; fi
  # call variants with haplotypecaller
  echo "Running HaplotypeCaller on scattered intervals for $interval"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo "$intervaldir/output.g.$scatter_interval.vcf.gz" >> /tmp/gvcf.list
    gatk --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap" HaplotypeCaller \
      -I $intervaldir/bqsr.$scatter_interval.bam \
      -R $reference/fasta/genome.fa \
      -O $intervaldir/output.g.$scatter_interval.vcf.gz \
      -ERC GVCF \
      -D $gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
      -L /tmp/interval_files_folder/$scatter_interval \
      --verbosity INFO \
      --tmp-dir $tmpdir >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33mHaplotypeCaller failed on $interval. Check $outputlog for additional info\033[0m"; echo 1 > /tmp/exit_status.txt; } &
  done
  # exit with 1 if interval failed
  wait
  exit_status=$(head -n1 /tmp/exit_status.txt)
  if [ $exit_status -eq 1 ]; then return 1; fi

  # merge and index vcfs
  gatk GatherVcfs -I /tmp/gvcf.list -O $intervaldir/output.g.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mGatherVcfs failed. Check $outputlog for additional info\033[0m"; return 1; }
  gatk IndexFeatureFile -I $intervaldir/output.g.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mIndexFeatureFile failed. Check $outputlog for additional info\033[0m"; return 1; }

  # single sample genotyping in parallel
  if [ -f /tmp/vcf.list ]; then rm /tmp/vcf.list; fi
  # genotype vcf
  echo "Genotyping scattered intervals for $interval"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo "$intervaldir/output.$scatter_interval.vcf.gz" >> /tmp/vcf.list
    gatk --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap" GenotypeGVCFs \
      -V $intervaldir/output.g.vcf.gz \
      -R $reference/fasta/genome.fa \
      -L /tmp/interval_files_folder/$scatter_interval \
      -O $intervaldir/output.$scatter_interval.vcf.gz >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33mGenotypeGVCFs failed on $interval. Check $outputlog for additional info\033[0m"; echo 1 > /tmp/exit_status.txt; } &
  done
  # exit with 1 if interval failed
  wait
  exit_status=$(head -n1 /tmp/exit_status.txt)
  if [ $exit_status -eq 1 ]; then return 1; fi

  # gather genotyped vcfs
  echo "Gathering genotyped intervals for $interval"
  gatk GatherVcfs -I /tmp/vcf.list -O $intervaldir/output.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mGatherVcfs failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }
  gatk IndexFeatureFile -I $intervaldir/output.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mIndexFeatureFile failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }
}

#####################################################################################################
#####################################################################################################
# RNA specific workflow
function interval_rna_germline_workflow {
  # workflow from here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
  echo "Starting RNA germline workflow for $interval"

  # if this step fails (and you are running docker in wsl) make sure there is enough disk space where docker is installed (usually C:/) and in $tmpdir
  # o/w you may see a cryptic jvm sigsegv samtools htslib error that results from generation of temp files in C:/users/$USER/AppData/Local/Docker/wsl
  # this may occur when the number of threads is increased (resulting in the generation of more temp files)
  if [ -f /tmp/cigarbam.list ]; then rm /tmp/cigarbam.list; fi
  echo "Running SplitNCigarReads on scattered intervals"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo "$intervaldir/cigar_marked_duplicates.$scatter_interval.bam" >> /tmp/cigarbam.list
    gatk --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap" SplitNCigarReads \
      -I $1 \
      -R $reference/fasta/genome.fa \
      --tmp-dir $tmpdir \
      -L /tmp/interval_files_folder/$scatter_interval \
      -O $intervaldir/cigar_marked_duplicates.$scatter_interval.bam >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33mSplitNCigarReads failed on $interval. Check $outputlog for additional info\033[0m"; echo 1 > /tmp/exit_status.txt; } &
  done
  # exit with 1 if interval failed
  wait
  exit_status=$(head -n1 /tmp/exit_status.txt)
  if [ $exit_status -eq 1 ]; then return 1; fi

  # gather the split cigar bam files and index
  echo "Gathering scattered bams for $interval"
  gatk GatherBamFiles \
    -I /tmp/cigarbam.list \
    -O $intervaldir/cigar_marked_duplicates.bam >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mGatherBamFiles failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }
  echo "Sorting gathered bam file"
  samtools sort -@ $threads $intervaldir/cigar_marked_duplicates.bam -o $intervaldir/sorted.cigar_marked_duplicates.bam >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33msamtools sort failed on $interval. Check $outputlog file for additional info\033[0m"; return 1; }
  samtools index -@ $threads $intervaldir/sorted.cigar_marked_duplicates.bam >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33msamtools index failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }

  # run gatk short variant pipeline using cigar split bam
  gatk_germline_short_variant_scatter_gather $intervaldir/sorted.cigar_marked_duplicates.bam

  # perform hard filtering using the qual-by-depth QD score and window for snp clustering
  # vqsr and cnnscorevariants is not recommended for rna-based genotyping
  # filter thresholds taken from https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl
  if [ -f /tmp/fvcf.list ]; then rm /tmp/fvcf.list; fi
  echo "Applying hard filters"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo $intervaldir/$library_id.$scatter_interval.vcf.gz >> /tmp/fvcf.list
    gatk VariantFiltration --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap" \
      --V $intervaldir/output.vcf.gz \
      --R $reference/fasta/genome.fa \
      -L /tmp/interval_files_folder/$scatter_interval \
      --verbosity ERROR \
      --window 35 \
      --cluster 3 \
      --filter-name "FS" \
      --filter "FS > 30.0" \
      --filter-name "QD" \
      --filter "QD < 2.0" \
      -O $intervaldir/$library_id.$scatter_interval.vcf.gz >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33mVariantFiltration failed on $interval. Check $outputlog for additional info\033[0m"; echo 1 > /tmp/exit_status.txt; } &
  done
  # exit with 1 if interval failed
  wait
  exit_status=$(head -n1 /tmp/exit_status.txt)
  if [ $exit_status -eq 1 ]; then return 1; fi

  # merge and index filtered vcfs
  echo "Gathering filtered intervals for $interval"
  gatk GatherVcfs -I /tmp/fvcf.list -O $workdir/genotype.$interval.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mGatherVcfs failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }
  gatk IndexFeatureFile -I $workdir/genotype.$interval.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mIndexFeatureFile failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }
} 

#####################################################################################################
#####################################################################################################
# ATAC specific workflow
function interval_atac_germline_workflow {
  echo "Starting ATAC germline workflow for $interval"
  # workflow from here: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
  # for single sample calling exclude joint-call cohort step
  # run variant pipeline using cellranger input bam
  gatk_germline_short_variant_scatter_gather $1

  # score the variants prior to filtering
  # cnnscorevariants may have better performance than vqsr for single sample genotyping. VQSR recommends >30 exomes
  # see discussion: https://gatk.broadinstitute.org/hc/en-us/community/posts/360056186812-Using-VQSR-for-small-scale-experiments
  if [ -f /tmp/cnnvcf.list ]; then rm /tmp/cnnvcf.list; fi
  echo "Running CNNScoreVariants on $interval"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo "$intervaldir/annotated.$scatter_interval.vcf.gz" >> /tmp/cnnvcf.list
    gatk CNNScoreVariants --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap" \
      -V $intervaldir/output.vcf.gz \
      -R $reference/fasta/genome.fa \
      -L /tmp/interval_files_folder/$scatter_interval \
      -O $intervaldir/annotated.$scatter_interval.vcf.gz >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33mCNNScoreVariants failed on $interval. Check $outputlog for additional info\033[0m"; echo 1 > /tmp/exit_status.txt; } &
  done
  # exit with 1 if interval failed
  wait
  exit_status=$(head -n1 /tmp/exit_status.txt)
  if [ $exit_status -eq 1 ]; then return 1; fi

  # gather cnn vcfs
  echo "Gathering CNNScoreVariants vcfs"
  gatk GatherVcfs -I /tmp/cnnvcf.list -O $intervaldir/annotated.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mGatherVcfs failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }
  gatk IndexFeatureFile -I $intervaldir/annotated.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mIndexFeatureFile failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }

  # filter variants with default tranches from gatk
  echo "Filtering gathered vcf using CNNScoreVariants tranches"
  gatk FilterVariantTranches --java-options "-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap" \
    -V $intervaldir/annotated.vcf.gz \
    --resource $gatk_bundle/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
    --resource $gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --info-key CNN_1D \
    --snp-tranche 99.95 \
    --indel-tranche 99.4 \
    -O $workdir/genotype.$interval.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mFilterVariantTranches failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }

  # index the interval vcf
  gatk IndexFeatureFile -I $workdir/genotype.$interval.vcf.gz >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mIndexFeatureFile failed on $interval. Check $outputlog for additional info\033[0m"; return 1; }
}

# check input files
###########################################################
###########################################################
# check if SCRATCH1 variable points to a directory
if [ ! -d $SCRATCH1 ]; then echo "$SCRATCH1 variable does not point to scratch directory"; exit 1; fi

# check for input bam file
if [ ! -f $inputbam ]; then echo "Input bam file not found"; exit 1; fi
if [ ! -f $inputbam.bai ]; then echo "Input bam index not found. Run samtools index on input bam"; exit 1; fi

# check for ref
if [ ! -f $reference/fasta/genome.fa ]; then
  echo "Reference genome.fa not found in $reference/fasta directory"
  exit 1
fi

# check for rna and atac gatk resource files for haplotype caller
# TODO: check for indexes and change names of files
if [ ! -f $gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz ] \
  || [ ! -f $gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz ] \
  || [ ! -f $gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz ] \
  || [ ! -f $gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ]; then
  echo "One or more GATK resource files not found in gatk_bundle directory:"
  echo "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz"
  echo "resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
  echo "resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  echo "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  exit 1
fi

# check for rna gtf file used to generate calling intervals
# TODO: make compatible with genes.gtf.gz file as seen in cellranger-arc reference
if [ $modality = rna ] && [ ! -f $reference/genes/genes.gtf ]; then
  echo "genes.gtf file not found in $reference/genes directory"
  exit 1
fi

# check for wgs calling regions and dna variant gatk resources
# TODO: check for indexes
if [ ! -f $gatk_bundle/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list ] \
  || [ ! -f $gatk_bundle/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz ] \
  || [ ! -f $gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ]; then
  if [ $modality = "atac" ]; then
  echo "One or more GATK resource files not found in gatk_bundle directory:"
  echo "resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
  echo "resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
  echo "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  exit 1
  fi
fi

##################################################################
##################################################################
# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# set exit status temporary file for monitoring return values in parallel processes
echo $exit_status > /tmp/exit_status.txt

# activate gatk conda environ
echo "Activating GATK conda environment"
source activate gatk

# specify a work directory
workdir=$SCRATCH1/gatk_genotype/$modality/$library_id
# create output directory
mkdir -p $outputdir
rm -rf $workdir; mkdir -p $workdir 2> /dev/null

# stream GATK output to terminal o/w capture in log file
if [ $verbose = "true" ]; then
  outputlog=/dev/stdout
elif [ $verbose = "false" ]; then
  outputlog=$workdir/log.out
fi

# prepare a fasta dict file using the cellranger ref if not already present
if [ ! -f $reference/fasta/genome.dict ] ; then
  gatk CreateSequenceDictionary -R $reference/fasta/genome.fa
fi

# create bed from gatk hg38 wgs calling regions for atac 
# or from gtf file used to create rna reference. Change gtf to 0-based coords by subtracting 1 from start in column 4
echo "Generating calling intervals bed file"
if [ $modality = "atac" ]; then
  grep -v @ $gatk_bundle/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list |pv| cut -f1-3 > /tmp/calling_intervals.bed
elif [ $modality = "rna" ]; then
# TODO: make compatible with genes.gtf.gz file as seen in cellranger-arc reference
  cat $reference/genes/genes.gtf| grep -v "#" |pv| awk -F'\t' 'BEGIN { OFS="\t" } $3=="exon" {print $1,$4-1,$5}' > /tmp/calling_intervals.bed
fi

# specify a temporary file directory for SplitNCigarReads and HaplotypeCaller
tmpdir=$workdir/Temp
mkdir -p $tmpdir

# create scatter gather intervals to break up alternating ACGT and N blocks
rm -rf /tmp/scatter_by_Ns.interval_list 2>> ${outputlog}
echo "Scattering intervals across reference"
gatk ScatterIntervalsByNs \
  -R $reference/fasta/genome.fa \
  -O /tmp/scatter_by_Ns.interval_list 

# limit to specified interval
if [ $interval ]; then
  echo "Selected chromosome is $interval"
  intervals=($interval)
else
  echo "Genotyping all intervals"
  intervals=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
fi

# genotype each chromosome in series
for interval in ${intervals[@]}; do
  intervaldir=$workdir/$interval
  mkdir $intervaldir

  # filter calling intervals by selected interval
  awk -v var=$interval -F'\t' 'BEGIN { OFS="\t" } $1==var {print $1,$2,$3}' /tmp/calling_intervals.bed > /tmp/calling_intervals_sel.bed

  # intersect scatter gather intervals with calling bed and divide across available no. threads
  rm -rf /tmp/interval_files_folder/ 2>> ${outputlog}
  gatk SplitIntervals \
    -R $reference/fasta/genome.fa \
    -O /tmp/interval_files_folder/ \
    --scatter-count $threads \
    --interval-set-rule INTERSECTION \
    -L /tmp/calling_intervals_sel.bed \
    -L /tmp/scatter_by_Ns.interval_list >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mSplitIntervals failed on $interval. Check $workdir/log.out for additional info\033[0m"; echo 1 > /tmp/exit_status.txt; exit 1; }

  # create array of scatter gather intervals  
  scatter_intervals=$(ls /tmp/interval_files_folder/)

  attempts=1
  if [ $modality = "rna" ]; then
    interval_rna_germline_workflow $inputbam || { echo -e "\033[0;31mRNA germline workflow failed for $interval on first attempt\033[0m"; attempts=$((attempts+1)); }
    if [ $attempts -eq 2 ]; then
      echo -e "\e[92mRestarting RNA germline workflow for $interval\033[0m"
      rm -rf $intervaldir; mkdir $intervaldir
      interval_rna_germline_workflow $inputbam || { echo -e "\033[0;31mRNA germline workflow failed for $interval on second attempt\033[0m"; echo 1 > /tmp/exit_status.txt; exit 1; }
    fi
  elif [ $modality = "atac" ]; then
    interval_atac_germline_workflow $inputbam || { echo -e "\033[0;31mATAC germline workflow failed for $interval on first attempt\033[0m"; attempts=$((attempts+1)); }
    if [ $attempts -eq 2 ]; then
      echo -e "\e[92mRestarting RNA germline workflow for $interval\033[0m"
      rm -rf $intervaldir; mkdir $intervaldir
      interval_atac_germline_workflow $inputbam || { echo -e "\033[0;31mATAC germline workflow failed for $interval on second attempt\033[0m"; echo 1 > /tmp/exit_status.txt; exit 1; }
    fi
  fi

  echo -e "\e[92mWriting genotyped vcf for $interval to $workdir \033[0m"

done | pv -t

# gather the genotyped vcf intervals
exit_status=$(head -n1 /tmp/exit_status.txt)
if [ $exit_status -eq 0 ]; then
  echo -e "\e[0;92mSaving $outputvcf to $outputdir \033[0m"
  ls -1 $workdir/genotype.chr*.vcf.gz > /tmp/final_vcf.list
  gatk MergeVcfs -I /tmp/final_vcf.list -O $outputdir/$outputvcf >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mGatherVcfs failed. Check $workdir/log.out for additional info\033[0m"; exit 1; }
    
  # save analysis-ready bam to outputdir  
  if [ $outputbam ]; then
    echo -e "\e[0;92mSaving $outputbam to $outputdir \033[0m"
    if [ -f /tmp/bam.list ]; then rm /tmp/bam.list; fi
    ls -1 $workdir/sorted.bqsr.chr*.bam > /tmp/bam.list
    gatk GatherBamFiles \
      -I /tmp/bam.list \
      -O $workdir/$library_id.$modality.bqsr.bam >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33mGatherBamFiles failed. Check $outputlog for additional info\033[0m"; return 1; }
    echo "Sorting gathered bam file"
    samtools sort -@ $threads $workdir/$library_id.$modality.bqsr.bam -o $outputdir/$outputbam >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33msamtools sort failed. Check $outputlog file for additional info\033[0m"; return 1; }
    samtools index -@ $threads $outputdir/$outputbam >> ${outputlog} 2>&1 \
      || { echo -e "\033[0;33msamtools index failed. Check $outputlog for additional info\033[0m"; return 1; }
  fi

  #cleanup
  rm -rf $workdir

  format_time() {
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "%02d:%02d:%02d\n" $h $m $s
   }

  echo -e "\033[0;40mGenotype completed in $(format_time $SECONDS)\033[0m"
fi
