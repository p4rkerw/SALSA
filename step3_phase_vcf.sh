#!/bin/bash

# Set some default values:
interval=""
threads=1
hcphase=false
snvonly=false
snvindel=false
reproduce=false

function usage {
        echo "Usage: $(basename $0) [-nvdolpsitrh]" 2>&1
        echo '   -n  | --library_id         STR   library_id: eg. [Control_1]'
        echo '   -v  | --inputvcf           STR   path/to/input.vcf.gz eg. [vcfdir/joint_genotype/Control_1.pass.joint.vcf.gz]'
        echo '   -d  | --outputdir          STR   output directory name eg. [vcfdir/phasing]'
        echo '   -o  | --outputvcf          STR   name of output vcf eg. [Control_1.pass.joint.phase.vcf.gz]'
        echo '   -l  | --interval           STR   optional: phase a single chromosome eg. [chr22]'
        echo '   -p  | --hcphase                  optional: recover haplotypecaller physical phasing variants that are not in shapeit reference. Default=[false]'
        echo '   -s  | --snvonly            STR   use the biallelic_SNV reference for phasing'
        echo '   -i  | --snvindel           STR   use the biallelic_SNV_and_INDEL reference for phasing'
        echo '   -r  | --reproduce                optional: run shapeit with a single thread for reproducibility. Default=[false]'
        echo '   -t  | --threads            INT   number of threads. Default=[1]'
        echo '   -h  | --help                     show usage'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

PARSED_ARGUMENTS=$(getopt -a -n step2_merge_genh.sh -o n:v:d:o:l:psirt:h --long library_id:,inputvcf:,outputdir:,outputvcf:,interval:,hcphase,snvonly,\
snvindel,reproduce,threads:,help -- "$@")

echo "PARSED_ARGUMENTS are $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -n | --library_id)        library_id=$2        ; shift 2 ;;
    -v | --inputvcf)          inputvcf=$2          ; shift 2 ;;
    -d | --outputdir)         outputdir=$2         ; shift 2 ;;
    -o | --outputvcf)         outputvcf=$2         ; shift 2 ;;
    -l | --interval)          interval=$2          ; shift 2 ;;
    -p | --hcphase)           hcphase=true         ; shift 1 ;;
    -s | --snvonly)           snvonly=true         ; shift 1 ;;
    -i | --snvindel)          snvindel=true        ; shift 1 ;;
    -r | --reproduce)         reproduce=true       ; shift 1 ;;
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
echo "interval                  : $interval"
echo "hcphase                   : $hcphase"
echo "snvonly                   : $snvonly"
echo "snvindel                  : $snvindel"
echo "reproduce                 : $reproduce"
echo "threads                   : $threads"
echo "Parameters remaining are  : $@"

# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# workflow for phasing vcf genotypes
mkdir -p $outputdir

workdir=$SCRATCH1/phasing/$library_id
rm -rf $workdir; mkdir -p $workdir 2> /dev/null

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

if [ ! -f $inputvcf.tbi ]; then
bcftools index --threads $threads --tbi $inputvcf
fi

# use a reference with both snv and indels for phasing
if [ $snvindel == "true" ] && [ $reproduce == "false" ]; then
  rm $workdir/vcf.list 2> /dev/null
  for interval in ${intervals[@]}; do
  echo "$workdir/phased.$interval.vcf.gz" >> $workdir/vcf.list
  shapeit4.2 \
  --thread $threads \
  --input $inputvcf \
  --map $workdir/phasing/shapeit4/maps/$interval.b38.gmap.gz \
  --region $interval \
  --seed 123456 \
  --reference phasing/ALL.$interval.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
  --sequencing \
  --output $workdir/phased.$interval.vcf.gz
  done
wait
fi

# reproduce flag set to true
if [ $snvindel == "true" ] && [ $reproduce == "true" ]; then
  rm $workdir/vcf.list 2> /dev/null
  for interval in ${intervals[@]}; do
  echo "$workdir/phased.$interval.vcf.gz" >> $workdir/vcf.list
  shapeit4.2 \
  --input $inputvcf \
  --map $workdir/phasing/shapeit4/maps/$interval.b38.gmap.gz \
  --region $interval \
  --seed 123456 \
  --reference phasing/ALL.$interval.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
  --sequencing \
  --output $workdir/phased.$interval.vcf.gz &
  done
wait
fi

# use a reference with only snv for phasing
if [ $snvonly = "true" ] && [ $reproduce == "false" ]; then
  rm $workdir/vcf.list 2> /dev/null
  for interval in ${intervals[@]}; do
  echo "$workdir/phased.$interval.vcf.gz" >> $workdir/vcf.list
  shapeit4.2 \
  --thread $threads \
  --input $inputvcf\
  --map $workdir/phasing/shapeit4/maps/$interval.b38.gmap.gz \
  --region $interval \
  --seed 123456 \
  --reference phasing/ALL.$interval.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz \
  --sequencing \
  --output $workdir/phased.$interval.vcf.gz
  done
  wait
fi

if [ $snvonly = "true" ] && [ $reproduce == "true" ]; then
  rm $workdir/vcf.list 2> /dev/null
  for interval in ${intervals[@]}; do
  echo "$workdir/phased.$interval.vcf.gz" >> $workdir/vcf.list
  shapeit4.2 \
  --input $inputvcf \
  --map $workdir/phasing/shapeit4/maps/$interval.b38.gmap.gz \
  --region $interval \
  --seed 123456 \
  --reference phasing/ALL.$interval.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz \
  --sequencing \
  --output $workdir/phased.$interval.vcf.gz &
  done
  wait
fi

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

  cp $workdir/anno_output.vcf.gz $outputdir/$outputvcf
  bcftools index --threads $threads --tbi $outputdir/$outputvcf
fi

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

echo -e "\033[35;40mPhasing completed in $(format_time $SECONDS)\033[0m"

# clean up
rm -rf $workdir