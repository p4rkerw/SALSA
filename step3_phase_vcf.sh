#!/bin/bash

# Set some default values:
interval=""
threads=1
hcphase=false
snvonly=false
snvindel=false
reproduce=false
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

Usage: step3_phase_vcf.sh [-nfdorlpsirvth]
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
  -V  | --verbose                  optional: stream shapeit4 output to terminal. Default=[false]
  -t  | --threads            INT   number of threads. Default=[1]
  -h  | --help                     show usage

EOF
exit 1
}

if [[ ${#} -eq 0 ]]; then usage; fi

PARSED_ARGUMENTS=$(getopt -a -n step3_phase_vcf.sh \
-o n:v:d:o:r:l:psirVt:h \
--long library_id:,inputvcf:,outputdir:,outputvcf:,phasingref:,interval:,hcphase,snvonly,snvindel,reproduce,verbose,threads:,help -- "$@")

echo "PARSED_ARGUMENTS are $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -n | --library_id)        library_id=$2        ; shift 2 ;;
    -v | --inputvcf)          inputvcf=$2          ; shift 2 ;;
    -d | --outputdir)         outputdir=$2         ; shift 2 ;;
    -o | --outputvcf)         outputvcf=$2         ; shift 2 ;;
    -r | --phasingref)        phasingref=$2        ; shift 2 ;;
    -l | --interval)          interval=$2          ; shift 2 ;;
    -p | --hcphase)           hcphase=true         ; shift 1 ;;
    -s | --snvonly)           snvonly=true         ; shift 1 ;;
    -i | --snvindel)          snvindel=true        ; shift 1 ;;
    -r | --reproduce)         reproduce=true       ; shift 1 ;;
    -V | --verbose)           verbose=true         ; shift 1 ;;
    -t | --threads)           threads=$2           ; shift 2 ;;
    -h | --help)              usage ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

echo "library_id                : $library_id"
echo "inputvcf                  : $inputvcf"
echo "outputdir                 : $outputdir"
echo "outputvcf                 : $outputvcf"
echo "phasingref                : $phasingref"
echo "interval                  : $interval"
echo "hcphase                   : $hcphase"
echo "snvonly                   : $snvonly"
echo "snvindel                  : $snvindel"
echo "reproduce                 : $reproduce"
echo "verbose                   : $verbose"
echo "threads                   : $threads"
echo "Parameters remaining are  : $@"

# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# check for inputs
if [ ! -f $inputvcf ]; then echo "Input vcf not found"; exit 1; fi

# TODO: check for reference phasing vcfs by interval if selected

# stream output to terminal o/w capture in log file
if [ $verbose = "true" ]; then
  outputlog=/dev/stdout
elif [ $verbose = "false" ]; then
  outputlog=$workdir/log.out
fi

# remove filtered variants
echo "Removing filtered variants from input vcf"
bcftools view -Oz -f PASS $inputvcf > /tmp/filter.vcf.gz
bcftools index --tbi --threads $threads /tmp/filter.vcf.gz

# workflow for phasing vcf genotypes
mkdir -p $outputdir

workdir=$SCRATCH1/phasing/$library_id
rm -rf $workdir; mkdir -p $workdir

# unpack individual chromosome maps
mkdir -p $workdir/phasing/shapeit4/maps
tar -xvzf /opt/shapeit4/maps/genetic_maps.b38.tar.gz -C $workdir/phasing/shapeit4/maps

# filter input vcf for selected intervals
if [ $interval ]; then
  echo "Selected interval is $interval"
  intervals=($interval)
else
  intervals=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
fi

# assign base shapeit args
shapeit_args=(--input /tmp/filter.vcf.gz \
  --seed 123456 \
  --sequencing)

# set multithreading if reproduce flag is false
if [ $reproduce = "false" ]; then
  shapeit_args+=(--thread ${threads})
fi

# remove vcf list from any prior phasing runs
rm $workdir/vcf.list 2> /dev/null

# set reference for SNV only or SNV_INDEL
if [ $snvindel = "true" ]; then
  for interval in ${intervals[@]}; do
  echo "$workdir/phased.$interval.vcf.gz" >> $workdir/vcf.list
  shapeit4.2 ${shapeit_args[@]} \
  --map $workdir/phasing/shapeit4/maps/${interval}.b38.gmap.gz \
  --region ${interval} \
  --reference $phasingref/ALL.${interval}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
  --output $workdir/phased.${interval}.vcf.gz >> ${outputlog}
  done
elif [ $snvonly = "true" ]; then
  for interval in ${intervals[@]}; do
  echo "$workdir/phased.$interval.vcf.gz" >> $workdir/vcf.list
  shapeit4.2 ${shapeit_args[@]} \
  --map $workdir/phasing/shapeit4/maps/${interval}.b38.gmap.gz \
  --region ${interval} \
  --reference $phasingref/ALL.${interval}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz \
  --output $workdir/phased.${interval}.vcf.gz >> ${outputlog}
  done
fi

# # reproduce flag set to false run in multithreaded series
# if [ $reproduce = "false" ]; then
#   pids=()
#   for interval in ${intervals[@]}; do
#     echo "$workdir/phased.$interval.vcf.gz" >> $workdir/vcf.list
#     shapeit4.2 ${shapeit_args[@]} >> ${outputlog}
#     pids+=($!)
#   done
#   # check exit status for each interva
#   for pid in ${pids[@]}; do
#     if ! wait $pid; then { echo "shapeit phasing failed check $outputlog"; exit_status=1; exit 1; };  fi
#   done
# fi

# # reproduce flag set to true run in parallel with single thread per interval
# if [ $reproduce = "true" ]; then
#   pids=()
#   for interval in ${intervals[@]}; do
#     echo "$workdir/phased.$interval.vcf.gz" >> $workdir/vcf.list
#     shapeit4.2 ${shapeit_args[@]} >> ${outputlog} &
#     pids+=($!)
#     while (( $(jobs |wc -l) >= (( ${threads} + 1 )) )); do
#       sleep 0.1
#     done
#   done
#   # check exit status for each interval
#   for pid in ${pids[@]}; do
#     if ! wait $pid; then { echo "shapeit phasing failed check $outputlog"; exit_status=1; exit 1; };  fi
#   done
# fi

# concatenate phased vcf files from all intervals
bcftools concat --threads $threads -Oz -f $workdir/vcf.list -o $outputdir/$outputvcf

# index the aggregated vcf
bcftools index --threads $threads --tbi $outputdir/$outputvcf

# ############################################
if [ $hcphase = "true" ]; then
  echo "Recovering physically phased variants that are not in shapeit reference"
  # # annotate shapeit phased variants with original vcf FORMAT field PID annotations from haplotypecaller
  bcftools annotate --threads $threads -Oz -c FORMAT,INFO -a $inputvcf $outputdir/$outputvcf > $workdir/shapeit_annotated.vcf.gz
  bcftools index --threads $threads $workdir/shapeit_annotated.vcf.gz

  # # extract PID tags associated with phased shapeit vcf
  bcftools filter -Oz -i 'FMT/PID!=""' $workdir/shapeit_annotated.vcf.gz | bcftools query -f'[%PID\n]'  > $workdir/shapeit_pid.txt

  # # extract all variants in input vcf that share a PID tag with a shapeit phased variant
  (bcftools view -h $inputvcf; bcftools view -H $inputvcf | grep -F -f $workdir/shapeit_pid.txt) | bcftools view -Oz - > $workdir/pid_variants_to_keep.vcf.gz
  bcftools index --threads $threads $workdir/pid_variants_to_keep.vcf.gz

  # # merge the variants that share a PID tag with a shapeit variant with the shapeit phased vcf
  bcftools concat --threads $threads -Oz --allow-overlaps --rm-dups "none" $workdir/shapeit_annotated.vcf.gz $workdir/pid_variants_to_keep.vcf.gz > $workdir/output.vcf.gz
  bcftools index --threads $threads $workdir/output.vcf.gz

  # # recover the genotype annotations in FORMAT field of shapeit phased vcf
  bcftools annotate --threads $threads -Oz -c FORMAT -a $outputdir/$outputvcf $workdir/output.vcf.gz > $workdir/anno_output.vcf.gz
  bcftools index --threads $threads $workdir/anno_output.vcf.gz

  # # # check output 
  bcftools query -f'[%CHROM,%POS,%REF,%ALT,%GT,%FILTER\n]' $workdir/anno_output.vcf.gz |head -n10
  bcftools query -f'[%GT\n]' $workdir/anno_output.vcf.gz |sort|uniq

  # # # select only phased variants and count records
  echo "Total number of variants in $inputvcf"
  bcftools view $inputvcf|bcftools stats --threads $threads -|grep record
  echo "Number of variants phased by shapeit"
  bcftools view $outputdir/$outputvcf|bcftools stats --threads $threads -|grep record
  echo "Number of variants phased by shapeit or physically phased by haplotypecaller"
  bcftools view -p $workdir/anno_output.vcf.gz|bcftools stats --threads $threads -|grep record

  mv $workdir/anno_output.vcf.gz $outputdir/$outputvcf
  bcftools index -f --threads $threads --tbi $outputdir/$outputvcf
fi

if [ $exit_status -eq 0 ]; then
  format_time() {
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "%02d:%02d:%02d\n" $h $m $s
   }

  echo -e "\e[92mWriting $outputvcf to $outputdir\033[0m"

  echo -e "\033[35;40mPhasing completed in $(format_time $SECONDS)\033[0m"

  # clean up
  rm -rf $workdir
fi