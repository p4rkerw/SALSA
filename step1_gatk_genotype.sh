#!/bin/bash
# this script will generate a genotyped vcf using the input bam file
# for info to implement a scatter workflow: https://www.ibm.com/downloads/cas/ZJQD0QAL

# exit on error
set -e

# Set some default values:
interval=""
threads=1
verbose="false"

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

Usage: step1_gatk_genotype.sh [-inrgdomlt]
  -i  | --inputbam           STR   path/to/input.bam eg. [project/sample_1/outs/possorted*.bam]
  -n  | --library_id         STR   library_id: eg. [sample_1]
  -r  | --reference          STR   path/to/cellranger_ref eg. [reference/refdata-gex-GRCh38-2020-A]
  -g  | --gatk_bundle        STR   path/to/gatk_bundle eg. [reference/gatk]
  -d  | --outputdir          STR   output directory name eg. [project/rna_genotype]
  -o  | --outputvcf          STR   name of output vcf eg. [sample_1.rna.vcf.gz]
  -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
  -l  | --interval           STR   optional: genotype a single chromosome eg. [chr22]
  -v  | --verbose                  optional: stream GATK output to terminal. Default=[false]
  -t  | --threads            INT   number of threads. Default=[1]
  -h  | --help                     show usage

EOF
exit 1
}

if [[ ${#} -eq 0 ]]; then usage; fi

PARSED_ARGUMENTS=$(getopt -a -n step1_gatk_genotype.sh \
-o i:n:r:g:d:o:m:l:vt:h \
--long inputbam:,library_id:,reference:,gatk_bundle:,outputdir:,outputvcf:,modality:,interval:,verbose,threads:,help -- "$@")

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
    -v | --verbose)             verbose=true                    ; shift 1 ;;
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
  pids=()
  for scatter_interval in ${scatter_intervals[@]}; do
    gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" BaseRecalibrator \
    -I $1 \
    -R $reference/fasta/genome.fa \
    --known-sites $gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --known-sites $gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites $gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites $gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -L /tmp/interval_files_folder/$scatter_interval \
    -O $workdir/recal_data.$scatter_interval.table \
    --verbosity INFO >> ${verbosity} 2>&1 \
      || { echo "BaseRecalibrator failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; } &
    pids+=($!)
  done 
  # check exit status for each interval
  for pid in ${pids[@]}; do
    if ! wait $pid; then exit 1; fi
  done

  # apply base quality score recalibration
  echo "Running ApplyBQSR on scattered intervals for $interval"
  for scatter_interval in ${scatter_intervals[@]}; do
    gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" ApplyBQSR \
    -I $1 \
    -R $reference/fasta/genome.fa \
    --bqsr-recal-file $workdir/recal_data.$scatter_interval.table \
    -L /tmp/interval_files_folder/$scatter_interval \
    -O $workdir/bqsr.$scatter_interval.bam \
    --verbosity INFO >> ${verbosity} 2>&1 \
      || { echo "ApplyBQSR failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; } &
  done
  wait

  # make sure bam outputs are indexed
  for scatter_interval in ${scatter_intervals[@]}; do
    if [ ! -f $workdir/bqsr.$scatter_interval.bam.bai ]; then
      samtools index -@ $threads $workdir/bqsr.$scatter_interval.bam
    fi
  done

  # remove vcf list if session has been used multiple times
  if [ -f /tmp/gvcf.list ]; then rm /tmp/gvcf.list; fi
  # call variants with haplotypecaller
  echo "Running HaplotypeCaller on scattered intervals for $interval"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo "$workdir/output.g.$scatter_interval.vcf.gz" >> /tmp/gvcf.list
    gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=1" HaplotypeCaller \
    -I $workdir/bqsr.$scatter_interval.bam \
    -R $reference/fasta/genome.fa \
    -O $workdir/output.g.$scatter_interval.vcf.gz \
    -ERC GVCF \
    -D $gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -L /tmp/interval_files_folder/$scatter_interval \
    --smith-waterman JAVA \
    --verbosity INFO \
    --tmp-dir $tmpdir >> ${verbosity} 2>&1 \
      || { echo "HaplotypeCaller failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; } &
  done
  wait

  # merge and index vcfs
  gatk GatherVcfs -I /tmp/gvcf.list -O $workdir/output.g.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "GatherVcfs failed. Check $SCRATCH1/log.out for additional info"; exit 1; }
  gatk IndexFeatureFile -I $workdir/output.g.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "IndexFeatureFile failed. Check $SCRATCH1/log.out for additional info"; exit 1; }

  # single sample genotyping in parallel
  if [ -f /tmp/vcf.list ]; then rm /tmp/vcf.list; fi
  # genotype vcf
  echo "Genotyping scattered intervals for $interval"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo "$workdir/output.$scatter_interval.vcf.gz" >> /tmp/vcf.list
    gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=1" GenotypeGVCFs \
    -V $workdir/output.g.vcf.gz \
    -R $reference/fasta/genome.fa \
    -L /tmp/interval_files_folder/$scatter_interval \
    -O $workdir/output.$scatter_interval.vcf.gz >> ${verbosity} 2>&1 \
      || { echo "GenotypeGVCFs failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; } &
  done
  wait

  # gather genotyped vcfs
  echo "Gathering genotyped intervals for $interval"
  gatk GatherVcfs -I /tmp/vcf.list -O $workdir/output.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "GatherVcfs failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }
  gatk IndexFeatureFile -I $workdir/output.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "IndexFeatureFile failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }
}

#####################################################################################################
#####################################################################################################
# RNA specific workflow
function interval_rna_germline_workflow {
  # workflow from here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
  inputbam=$1
  echo "Starting RNA germline workflow for $interval"

  # if this step fails (and you are running docker in wsl) make sure there is enough disk space where docker is installed (usually C:/) and in $tmpdir
  # o/w you may see a cryptic jvm sigsegv samtools htslib error that results from generation of temp files in C:/users/$USER/AppData/Local/Docker/wsl
  # this may occur when the number of threads is increased (resulting in the generation of more temp files)
  if [ -f /tmp/cigarbam.list ]; then rm /tmp/cigarbam.list; fi
  echo "Running SplitNCigarReads on scattered intervals"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo "$workdir/cigar_marked_duplicates.$scatter_interval.bam" >> /tmp/cigarbam.list
    gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=1" SplitNCigarReads \
    -I $inputbam \
    -R $reference/fasta/genome.fa \
    --tmp-dir $tmpdir \
    -L /tmp/interval_files_folder/$scatter_interval \
    -O $workdir/cigar_marked_duplicates.$scatter_interval.bam >> ${verbosity} 2>&1 \
      || { echo "SplitNCigarReads failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; } &
  done
  wait

  # gather the split cigar bam files and index
  echo "Gathering scattered bams for $interval"
  gatk GatherBamFiles \
  -I /tmp/cigarbam.list \
  -O $workdir/cigar_marked_duplicates.bam >> ${verbosity} 2>&1 \
    || { echo "GatherBamFiles failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }
  echo "Sorting gathered bam file"
  samtools sort -@ $threads $workdir/cigar_marked_duplicates.bam -o $workdir/sorted.cigar_marked_duplicates.bam >> ${verbosity} 2>&1 \
    || { echo "samtools sort failed on $interval. Check $SCRATCH1/log.out file for additional info"; exit 1; }
  samtools index -@ $threads $workdir/sorted.cigar_marked_duplicates.bam >> ${verbosity} 2>&1 \
    || { echo "samtools index failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }

  # run gatk short variant pipeline using cigar split bam
  gatk_germline_short_variant_scatter_gather $workdir/sorted.cigar_marked_duplicates.bam

  # perform hard filtering using the qual-by-depth QD score and window for snp clustering
  # vqsr and cnnscorevariants is not recommended for rna-based genotyping
  # filter thresholds taken from https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl
  if [ -f /tmp/fvcf.list ]; then rm /tmp/fvcf.list; fi
  echo "Applying hard filters"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo $workdir/$library_id.$scatter_interval.vcf.gz >> /tmp/fvcf.list
    gatk VariantFiltration \
    --V $workdir/output.vcf.gz \
    --R $reference/fasta/genome.fa \
    -L /tmp/interval_files_folder/$scatter_interval \
    --verbosity ERROR \
    --window 35 \
    --cluster 3 \
    --filter-name "FS" \
    --filter "FS > 30.0" \
    --filter-name "QD" \
    --filter "QD < 2.0" \
    -O $workdir/$library_id.$scatter_interval.vcf.gz >> ${verbosity} 2>&1 \
      || { echo "VariantFiltration failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; } &
  done
  wait

  # merge and index filtered vcfs
  echo "Gathering filtered intervals for $interval"
  gatk GatherVcfs -I /tmp/fvcf.list -O $workdir/genotype.$interval.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "GatherVcfs failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }
  gatk IndexFeatureFile -I $workdir/genotype.$interval.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "IndexFeatureFile failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }
} 

#####################################################################################################
#####################################################################################################
# ATAC specific workflow
function interval_atac_germline_workflow {
  echo "Starting ATAC germline workflow for $interval"

  inputbam=$1

  # before running pipeline ensure interval is coordinate-sorted and @HD tag is properly formatted
  samtools view -h $inputbam $interval > $workdir/$interval.bam
  # samtools sort -T $workdir -@ $threads $workdir/$interval.bam > $workdir/sorted.$interval.bam
  samtools index -@ $threads $workdir/$interval.bam

  # workflow from here: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
  # for single sample calling exclude joint-call cohort step
  # run variant pipeline using cellranger input bam
  gatk_germline_short_variant_scatter_gather $workdir/$interval.bam

  # score the variants prior to filtering
  # cnnscorevariants may have better performance than vqsr for single sample genotyping. VQSR recommends >30 exomes
  # see discussion: https://gatk.broadinstitute.org/hc/en-us/community/posts/360056186812-Using-VQSR-for-small-scale-experiments
  if [ -f /tmp/cnnvcf.list ]; then rm /tmp/cnnvcf.list; fi
  echo "Runnning CNNScoreVariants on $interval"
  for scatter_interval in ${scatter_intervals[@]}; do
    echo "$workdir/annotated.$scatter_interval.vcf.gz" >> /tmp/cnnvcf.list
    gatk CNNScoreVariants \
    -V $workdir/output.vcf.gz \
    -R $reference/fasta/genome.fa \
    -L /tmp/interval_files_folder/$scatter_interval \
    -O $workdir/annotated.$scatter_interval.vcf.gz >> ${verbosity} 2>&1 \
      || { echo "CNNScoreVariants failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; } &
  done
  wait

  # gather cnn vcfs
  echo "Gathering CNNScoreVariants vcfs"
  gatk GatherVcfs -I /tmp/cnnvcf.list -O $workdir/annotated.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "GatherVcfs failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }
  gatk IndexFeatureFile -I $workdir/annotated.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "IndexFeatureFile failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }

  # filter variants with default tranches from gatk
  echo "Filtering gathered vcf using CNNScoreVariants tranches"
  gatk FilterVariantTranches \
  -V $workdir/annotated.vcf.gz \
  --resource $gatk_bundle/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
  --resource $gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --info-key CNN_1D \
  --snp-tranche 99.95 \
  --indel-tranche 99.4 \
  -O $workdir/genotype.$interval.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "FilterVariantTranches failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }

  # index the interval vcf
  gatk IndexFeatureFile -I $workdir/genotype.$interval.vcf.gz >> ${verbosity} 2>&1 \
    || { echo "IndexFeatureFile failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }
}

# check input files
###########################################################
###########################################################
# check for input bam file
if [ ! -f $inputbam ]; then echo "Input bam file not found"; exit 1; fi

# check for ref
if [ ! -f $reference/fasta/genome.fa ]; then
  echo "Reference genome.fa not found in $reference/fasta directory"
  exit 1
fi

# check for rna and atac gatk resource files for haplotype caller
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
if [ $modality = rna ] && [ ! -f $reference/genes/genes.gtf ]; then
  echo "genes.gtf file not found in $reference/genes directory"
  exit 1
fi

# check for wgs calling regions and dna variant gatk resources
if [ ! -f $gatk_bundle/resources_broad_hg38_v0_wgs_calling_regions.hg38 ] \
  || [ ! -f $gatk_bundle/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz ] \
  || [ ! -f $gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ]; then
  if [ $modality = "atac" ]; then
  echo "One or more GATK resource files not found in gatk_bundle directory:"
  echo "resources_broad_hg38_v0_wgs_calling_regions.hg38"
  echo "resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
  echo "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  exit 1
  fi
fi

##################################################################
##################################################################
# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# activate gatk conda environ
echo "Activating GATK conda environment"
source activate gatk

# stream GATK output to terminal o/w capture in log file
if [ $verbose == "true" ]; then
  verbosity=/dev/stdout
elif [ $verbose == "false" ]; then
  verbosity=$SCRATCH1/log.out
fi

# prepare a fasta dict file using the cellranger ref if not already present
if [ ! -f $reference/fasta/genome.dict ] ; then
   gatk CreateSequenceDictionary -R $reference/fasta/genome.fa
fi

# create bed from gatk hg38 wgs calling regions for atac 
# or from gtf file used to create rna reference. Change gtf to 0-based coords by subtracting 1 from start in column 4
echo "Generating calling intervals bed file"
if [ $modality == "atac" ]; then
  grep -v @ $gatk_bundle/resources_broad_hg38_v0_wgs_calling_regions.hg38 |pv| cut -f1-3 > /tmp/calling_intervals.bed
elif [ $modality == "rna" ]; then
  grep -v '#' $reference/genes/genes.gtf |pv| awk -F'\t' 'BEGIN { OFS="\t" } $3=="exon" {print $1,$4-1,$5}' > /tmp/calling_intervals.bed
fi

# create a output and work directories
mkdir -p $outputdir
workdir=$SCRATCH1/gatk_genotype/$modality/$library_id
rm -rf $workdir; mkdir -p $workdir 2>> ${verbosity}

# specify a temporary file directory for SplitNCigarReads and HaplotypeCaller
tmpdir=$workdir/Temp
mkdir -p $tmpdir

# create scatter gather intervals to break up alternating ACGT and N blocks
rm -rf /tmp/scatter_by_Ns.interval_list 2>> ${verbosity}
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

  # filter calling intervals by selected interval
  awk -v var=$interval -F'\t' 'BEGIN { OFS="\t" } $1==var {print $1,$2,$3}' /tmp/calling_intervals.bed > /tmp/calling_intervals_sel.bed

  # intersect scatter gather intervals with calling bed and divide across available no. threads
  rm -rf /tmp/interval_files_folder/ 2>> ${verbosity}
  gatk SplitIntervals \
  -R $reference/fasta/genome.fa \
  -O /tmp/interval_files_folder/ \
  --scatter-count $threads \
  --interval-set-rule INTERSECTION \
  -L /tmp/calling_intervals_sel.bed \
  -L /tmp/scatter_by_Ns.interval_list >> ${verbosity} 2>&1 \
    || { echo "SplitIntervals failed on $interval. Check $SCRATCH1/log.out for additional info"; exit 1; }

  # create array of scatter gather intervals  
  scatter_intervals=$(ls /tmp/interval_files_folder/)

  if [ $modality == "rna" ]; then
    interval_rna_germline_workflow $inputbam
  elif [ $modality == "atac" ]; then
    interval_atac_germline_workflow $inputbam
  fi

  echo -e "\e[92mWriting genotyped vcf for $interval to $workdir \033[0m"

done | pv -t


# gather the genotyped vcf intervals
echo "Saving $outputvcf to $outputdir"
ls -1 $workdir/genotype.chr*.vcf.gz > /tmp/final_vcf.list
gatk MergeVcfs -I /tmp/final_vcf.list -O $outputdir/$outputvcf >> ${verbosity} 2>&1 \
  || { echo "GatherVcfs failed. Check $SCRATCH1/log.out for additional info"; exit 1; }


#cleanup
rm -rf $workdir

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

echo -e "\033[35;40mGenotype completed in $(format_time $SECONDS)\033[0m"