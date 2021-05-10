#!/bin/bash

# Set some default values:
includefiltered=false
threads=1

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

Usage: step2_merge_geno.sh [-nabdit]
  -n  | --library_id         STR   library_id: eg. [sample_1]
  -a  | --vcfone             STR   path/to/first.vcf.gz eg. [project/atac_genotype/sample_1.atac.vcf.gz]. Prioritize annotations from first vcf if there is overlap.
  -b  | --vcftwo             STR   path/to/second.vcf.gz eg. [project/rna_genotype/sample_1.rna.vcf.gz]
  -d  | --outputdir          STR   output directory [project/joint_genotype]
  -o  | --outputvcf          STR   name of output vcf eg. [sample_1.pass.joint.vcf.gz]
  -i  | --includefiltered          optional: include filtered variants in output vcf. Default=[false]
  -t  | --threads            INT   number of threads. Default=[1]
  -h  | --help                     show usage

EOF
exit 1
}


if [[ ${#} -eq 0 ]]; then
   usage
fi

PARSED_ARGUMENTS=$(getopt -a -n step2_merge_geno.sh \
-o n:a:b:d:o:it:h \
--long library_id:,vcfone:,vcftwo:,outputdir:,outputvcf:,includefiltered,threads:,help -- "$@")

echo "PARSED_ARGUMENTS are $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -n | --library_id)            library_id=$2        ; shift 2 ;;
    -a | --vcfone)                vcfone=$2            ; shift 2 ;;
    -b | --vcftwo)                vcftwo=$2            ; shift 2 ;;
    -d | --outputdir)             outputdir=$2         ; shift 2 ;;
    -o | --outputvcf)             outputvcf=$2         ; shift 2 ;;
    -i | --includefiltered)       includefiltered=true ; shift 1 ;;
    -t | --threads)               threads=$2           ; shift 2 ;;
    -h | --help)                  usage ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

echo "library_id                : $library_id"
echo "vcfone                    : $vcfone"
echo "vcftwo                    : $vcftwo"
echo "outputdir                 : $outputdir"
echo "outputvcf                 : $outputvcf"
echo "includefiltered           : $includefiltered"
echo "threads                   : $threads"
echo "Parameters remaining are  : $@"

# check input vcfs
if [ ! -f $vcfone ] || [ ! -f $vcftwo ]; then { echo "Input vcfs not found"; exit 1; }; fi

# start script
echo "Merging vcf files:"
echo $vcfone 
echo $vcftwo

mkdir vcfdir/$outputdir 2> /dev/null

# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# include filtered variants in output if indicated o/w remove (default)
# prioritize format field annotations in vcfone
if [ $includefiltered = "true" ]; then
	# identify intersecting variants and retain annotation from first vcf if there is overlap
	# ie add variants unique to vcftwo to variants in vcfone
	# use bcftools concat/norm to join biallelic snps and indels into multiallelic records which are then filtered. retain single info/format field
	echo "Retaining filtered variants"
	bcftools isec -Oz --threads $threads -p /tmp/isec $vcfone $vcftwo
	bcftools concat --threads $threads --allow-overlaps --rm-dups "none" /tmp/isec/0000.vcf.gz /tmp/isec/0001.vcf.gz /tmp/isec/0002.vcf.gz \
	  | bcftools norm -m +both |bcftools view -Oz -m2 -M2 > $outputdir/$outputvcf
	bcftools index --threads $threads $outputdir/$outputvcf
elif [ $includefiltered = "false" ]; then
	echo "Removing filtered variants"
	# remove filtered variants
	bcftools isec -f PASS -Oz --threads $threads -p /tmp/isec $vcfone $vcftwo
	bcftools concat --threads $threads --allow-overlaps --rm-dups "none" /tmp/isec/0000.vcf.gz /tmp/isec/0001.vcf.gz /tmp/isec/0002.vcf.gz \
	  | bcftools norm -m +both |bcftools view -Oz -f PASS -m2 -M2 > $outputdir/$outputvcf
	bcftools index --threads $threads $outputdir/$outputvcf
fi

echo -e "\e[92mWriting $outputvcf to $outputdir\033[0m"
