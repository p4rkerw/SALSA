#!/bin/bash
# this script will generate a genotype vcf using the input bam file
# parallelization is currently implemented only for haplotypecaller
# for info to implement a scatter workflow: https://www.ibm.com/downloads/cas/ZJQD0QAL

# exit on error
set -e

# Set some default values:
interval=""
threads=1

function usage {
        echo "Usage: $(basename $0) [-indomlt]" 2>&1
        echo '   -i  | --bam                STR   path/to/input.bam eg. [rna_counts/Control_1/outs/possorted*.bam]'
        echo '   -n  | --library_id         STR   library_id: eg. [Control_1]'
        echo '   -d  | --outputdir          STR   output directory name eg. [vcfdir/rna_genotype]'
        echo '   -o  | --outputvcf          STR   name of output vcf eg. [Control_1.rna.vcf.gz]'
        echo '   -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]'
        echo '   -l  | --interval           STR   optional: genotype a single chromosome eg. [chr10]'
        echo '   -t  | --threads            INT   number of threads. Default=[1]'
        echo '   -h  | --help                     show usage'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

PARSED_ARGUMENTS=$(getopt -a -n step1_gatk_genotype.sh -o i:n:d:o:m:l:t:h --long bam:,library_id:,outputdir:,outputvcf:,modality:,interval:,threads:,help -- "$@")

echo "PARSED_ARGUMENTS are $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -i | --bam)                 inputbam=$2                     ; shift 2 ;;
    -n | --library_id)          library_id=$2                   ; shift 2 ;;
    -d | --outputdir)           outputdir=$2                    ; shift 2 ;;
    -o | --outputvcf)           outputvcf=$2                    ; shift 2 ;;
    -m | --modality)            modality=$2                     ; shift 2 ;;
    -l | --interval)            interval=$2                     ; shift 2 ;;
    -t | --threads)             threads=$2                      ; shift 2 ;;
    -h | --help)                usage ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

echo "inputbam                  : $inputbam"
echo "library_id                : $library_id"
echo "outputdir                 : $outputdir"
echo "outputvcf                 : $outputvcf"
echo "modality                  : $modality"
echo "interval                  : $interval"
echo "threads                   : $threads"
echo "Parameters remaining are  : $@"

# check input files
###########################################################
# check for input bam file
if [ ! -f $inputbam ]; then
  echo "Input bam file not found"
  exit 1
fi

# check for ref
if [ ! -f ${modality}_ref/fasta/genome.fa ]; then
  echo "Reference genome.fa not found in ${modality}_ref/fasta directory"
  exit 1
fi

# check for rna and atac gatk resource files for haplotype caller
if [ ! -f gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz ] ||\
   [ ! -f gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz ] ||\
   [ ! -f gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz ] ||\
   [ ! -f gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ]; then
  echo "gatk resource files not found in gatk_bundle directory:"
  echo "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz"
  echo "resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
  echo "resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  echo "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  exit 1
fi

# check for rna gtf file used to generate calling intervals
if [ $modality = rna ] && [ ! -f rna_ref/genes/genes.gtf ]; then
  echo "genes.gtf file not found in rna_ref/genes directory"
  exit 1
fi

# check for wgs calling regions and dna variant gatk resources
if [ ! -f gatk_bundle/resources_broad_hg38_v0_wgs_calling_regions.hg38 ] ||\
   [ ! -f gatk_bundle/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz ] ||\
   [ ! -f gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ]; then
   if [ $modality = "atac" ]; then
    echo "gatk resource files not found in gatk_bundle directory:"
    echo "resources_broad_hg38_v0_wgs_calling_regions.hg38"
    echo "resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
    echo "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    exit 1
   fi
fi

# check if vcfdir is properly mounted
if [ ! -d vcfdir ]; then
  echo "vcf output directory not found. Mount a docker volume to $HOME/vcfdir"
fi

##################################################################
# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# activate gatk conda environ
source activate gatk

# prepare a fasta dict file using the cellranger ref if not already present
if [ ! -f ${modality}_ref/fasta/genome.dict ] ; then
   gatk CreateSequenceDictionary -R ${modality}_ref/fasta/genome.fa
fi

# create an output directory
mkdir -p $outputdir

# create bed from gatk hg38 wgs calling regions for atac 
# or from gtf file used to create rna reference. Change gtf to 0-based coords by subtracting 1 from start in column 4
if [ $modality = atac ]; then
  grep -v @ gatk_bundle/resources_broad_hg38_v0_wgs_calling_regions.hg38 | cut -f1-3 > /tmp/calling_intervals.bed
elif [ $modality = rna ]; then
  grep -v '#' rna_ref/genes/genes.gtf | awk -F'\t' 'BEGIN { OFS="\t" } $3=="exon" {print $1,$4-1,$5}' > /tmp/calling_intervals.bed
fi

# set calling regions
# grep -v @ gatk_bundle/resources_broad_hg38_v0_wgs_calling_regions.hg38 | cut -f1-3 > /tmp/calling_intervals.bed

# limit to specified intervals for testing purposes if -l flag selected
if [ $interval ]; then
  echo "Selected interval is $interval"
  verbosity=INFO
  workdir=$SCRATCH1/gatk_genotype/${modality}/$library_id/$interval
  rm -rf $workdir; mkdir -p $workdir 2> /dev/null

  # filter calling intervals by selected interval
  awk -v var=$interval -F'\t' 'BEGIN { OFS="\t" } $1==var {print $1,$2,$3}' /tmp/calling_intervals.bed > /tmp/calling_intervals_sel.bed
  mv /tmp/calling_intervals_sel.bed /tmp/calling_intervals.bed

  # filter cellranger bam for selected interval
  # echo "Filtering cellranger bam for selected interval $interval"
  # samtools view -@ $threads -bS $inputbam $interval|pv|samtools sort -@ $threads > $workdir/$interval.filtered.bam 
  # inputbam=$workdir/$interval.filtered.bam
  # samtools index -@ $threads $inputbam
else
  verbosity=ERROR      
  workdir=$SCRATCH1/gatk_genotype/${modality}/$library_id/allcontigs
  rm -rf $workdir; mkdir -p $workdir 2> /dev/null
fi

# create a temporary file directory
tmpdir=$workdir/Temp
mkdir -p $tmpdir

# create scatter gather intervals to break up alternating ACGT and N blocks for haplotypecaller
rm -rf /tmp/scatter_by_Ns.interval_list 2> /dev/null
gatk ScatterIntervalsByNs \
-R ${modality}_ref/fasta/genome.fa \
-O /tmp/scatter_by_Ns.interval_list

# intersect scatter gather intervals with calling bed and divide across available no. threads
# if interval is specified limit intervals to interval_sel
rm -rf /tmp/interval_files_folder 2> /dev/null
gatk SplitIntervals \
-R ${modality}_ref/fasta/genome.fa \
-O /tmp/interval_files_folder \
--scatter-count $threads \
--interval-set-rule INTERSECTION \
-L /tmp/calling_intervals.bed \
-L /tmp/scatter_by_Ns.interval_list

# create array of scatter gather intervals  
intervals=$(ls /tmp/interval_files_folder)

function gatk_germline_short_variant_scatter_gather {
# input bam is first positional argument $1
# bundle files can be obtained from gatk resource bundle on google cloud
# split selected intervals across the number of threads
# generate base recalibration table
  for interval in ${intervals[@]}; do
    gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" BaseRecalibrator \
    -I $1 \
    -R ${modality}_ref/fasta/genome.fa \
    --known-sites gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --known-sites gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -L /tmp/interval_files_folder/$interval \
    -O $workdir/recal_data.$interval.table \
    --verbosity INFO &
  done
  wait

# # # apply base quality score recalibration
  for interval in ${intervals[@]}; do
    gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" ApplyBQSR \
    -I $1 \
    -R ${modality}_ref/fasta/genome.fa \
    --bqsr-recal-file $workdir/recal_data.$interval.table \
    -L /tmp/interval_files_folder/$interval \
    -O $workdir/bqsr.$interval.bam \
    --verbosity INFO &
  done
  wait

  # remove vcf list if session has been used multiple times
  if [ -f /tmp/gvcf.list ]; then
    rm /tmp/gvcf.list
  fi
  # call variants with haplotypecaller
  for interval in ${intervals[@]}; do
    echo "$workdir/output.g.$interval.vcf.gz" >> /tmp/gvcf.list
    gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=1" HaplotypeCaller \
    -I $workdir/bqsr.$interval.bam \
    -R ${modality}_ref/fasta/genome.fa \
    -O $workdir/output.g.$interval.vcf.gz \
    -ERC GVCF \
    -D gatk_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -L /tmp/interval_files_folder/$interval \
    --verbosity INFO \
    --tmp-dir $tmpdir &
  done
  wait

  # merge and index vcfs
  gatk GatherVcfs -I /tmp/gvcf.list -O $workdir/output.g.vcf.gz
  gatk IndexFeatureFile -I $workdir/output.g.vcf.gz

  # single sample genotyping in parallel
  if [ -f /tmp/vcf.list ]; then
    rm /tmp/vcf.list
  fi
  # genotype vcf
  for interval in ${intervals[@]}; do
    echo "$workdir/output.$interval.vcf.gz" >> /tmp/vcf.list
    gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=1" GenotypeGVCFs \
    -V $workdir/output.g.vcf.gz \
    -R ${modality}_ref/fasta/genome.fa \
    -L /tmp/interval_files_folder/$interval \
    -O $workdir/output.$interval.vcf.gz &
  done
  wait

  # gather genotyped vcfs
  gatk GatherVcfs -I /tmp/vcf.list -O $workdir/output.vcf.gz
  gatk IndexFeatureFile -I $workdir/output.vcf.gz
}

#####################################################################################################
#####################################################################################################
# RNA specific workflow
if [ $modality = rna ]; then
  # workflow from here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
  echo $inputbam

  # splitcigar for selected intervals o/w alt contigs and unmapped reads are also processed
  # if this step fails (and you are running docker in wsl) make sure there is enough disk space where docker is installed (usually C:/) and in $tmpdir
  # o/w you may see a cryptic jvm sigsegv samtools htslib error that results from generation of temp files in C:/users/$USER/AppData/Local/Docker/wsl
  # this may occur when the number of threads is increased (resulting in the generation of more temp files)
  if [ -f /tmp/cigarbam.list ]; then
    rm /tmp/cigarbam.list
  fi
  for interval in ${intervals[@]}; do
    echo "$workdir/cigar_marked_duplicates.$interval.bam" >> /tmp/cigarbam.list
    gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" SplitNCigarReads \
    -I $inputbam \
    -R ${modality}_ref/fasta/genome.fa \
    --tmp-dir $tmpdir \
    -L /tmp/interval_files_folder/$interval \
    -O $workdir/cigar_marked_duplicates.$interval.bam &
  done
  wait

  # gather the split cigar bam files and index
  gatk GatherBamFiles \
  -I /tmp/cigarbam.list \
  -O $workdir/cigar_marked_duplicates.bam
  samtools sort -@ $threads $workdir/cigar_marked_duplicates.bam -o $workdir/sorted.cigar_marked_duplicates.bam
  samtools index -@ $threads $workdir/sorted.cigar_marked_duplicates.bam

  # run gatk short variant pipeline using cigar split bam
  gatk_germline_short_variant_scatter_gather $workdir/sorted.cigar_marked_duplicates.bam

  # perform hard filtering using the qual-by-depth QD score and window for snp clustering
  # vqsr and cnnscorevariants is not recommended for rna-based genotyping
  # filter thresholds taken from https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl
  if [ -f /tmp/fvcf.list ]; then
    rm /tmp/fvcf.list
  fi
  for interval in ${intervals[@]}; do
    echo $workdir/$library_id.$interval.vcf.gz >> /tmp/fvcf.list
    gatk VariantFiltration \
    --V $workdir/output.vcf.gz \
    --R ${modality}_ref/fasta/genome.fa \
    -L /tmp/interval_files_folder/$interval \
    --verbosity ERROR \
    --window 35 \
    --cluster 3 \
    --filter-name "FS" \
    --filter "FS > 30.0" \
    --filter-name "QD" \
    --filter "QD < 2.0" \
    -O $workdir/$library_id.$interval.vcf.gz &
  done
  wait

  # merge and index filtered vcfs
  gatk GatherVcfs -I /tmp/fvcf.list -O $outputdir/$outputvcf
  gatk IndexFeatureFile -I $outputdir/$outputvcf
fi 

#####################################################################################################
#####################################################################################################
# ATAC specific workflow
if [ $modality = atac ]; then
  # workflow from here: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
  # for single sample calling exclude joint-call cohort step
  # run variant pipeline using cellranger input bam
  gatk_germline_short_variant_scatter_gather $inputbam

  # score the variants prior to filtering
  # activate the conda packages that are required for cnnscorevariants
  # cnnscorevariants may have better performance than vqsr for single sample genotyping. VQSR recommends >30 exomes
  # see discussion: https://gatk.broadinstitute.org/hc/en-us/community/posts/360056186812-Using-VQSR-for-small-scale-experiments
  if [ -f /tmp/cnnvcf.list ]; then
    rm /tmp/cnnvcf.list
  fi
  # genotype vcf
  for interval in ${intervals[@]}; do
    echo "$workdir/annotated.$interval.vcf.gz" >> /tmp/cnnvcf.list
    gatk CNNScoreVariants \
    -V $workdir/output.vcf.gz \
    -R ${modality}_ref/fasta/genome.fa \
    -L /tmp/interval_files_folder/$interval \
    -O $workdir/annotated.$interval.vcf.gz &
  done
  wait

  # gather cnn vcfs
  gatk GatherVcfs -I /tmp/cnnvcf.list -O $workdir/annotated.vcf.gz
  gatk IndexFeatureFile -I $workdir/annotated.vcf.gz

  # filter variants with default tranches from gatk
  gatk FilterVariantTranches \
  -V $workdir/annotated.vcf.gz \
  --resource gatk_bundle/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
  --resource gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --info-key CNN_1D \
  --snp-tranche 99.95 \
  --indel-tranche 99.4 \
  -O $outputdir/$outputvcf
fi  

#cleanup
rm -rf $workdir

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

echo -e "\033[35;40mGenotype completed in $(format_time $SECONDS)\033[0m"