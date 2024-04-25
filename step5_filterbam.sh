#!/bin/bash
# this script will take a cellranger bam and filter for reads with valid barcode tags to prepare for WASP

# default values
threads=1
interval=""
validate=false

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

Usage: step5_filterbam.sh [-nidolmbeth]
-n  | --library_id         STR   library_id: eg. [sample_1]
-i  | --inputbam           STR   path/to/input.bam eg. [project/sample_1/outs/possorted_genome_bam.bam]
-d  | --outputdir          STR   output directory eg. [project/wasp_rna]
-o  | --outputbam          STR   filtered output bam eg. [sample_1.bcfilter.bam]
-l  | --interval           STR   optional: filter a single chromosome eg. [chr22]
-m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]
-b  | --barcodes           STR   path/to/barcodes.csv with headers and three columns. First column is named "barcode"
                                 second column is group "orig.ident" and third column is "celltype" eg. [project/barcodes/rna_barcodes.csv]
-e  | --validate                 validate the barcode-filtered bam file. Default=[false]
-t  | --threads            INT   number of threads. Default=[1]
-h  | --help                     show usage

EOF
exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

PARSED_ARGUMENTS=$(getopt -a -n step5_filterbam.sh \
-o n:i:d:o:l:m:b:e:t:h \
--long library_id:,inputbam:,outputdir:,outputbam:,interval:,modality:,barcodes:,validate,threads:,help -- "$@")

echo "PARSED_ARGUMENTS are $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -n | --library_id)        library_id=$2        ; shift 2 ;;
    -i | --inputbam)          inputbam=$2          ; shift 2 ;;
    -d | --outputdir)         outputdir=$2         ; shift 2 ;;
    -o | --outputbam)         outputbam=$2         ; shift 2 ;;
    -l | --interval)          interval=$2          ; shift 2 ;;
    -m | --modality)          modality=$2          ; shift 2 ;;
    -b | --barcodes)          barcodes=$2          ; shift 2 ;;
    -e | --validate)          validate=true        ; shift 1 ;;
    -t | --threads)           threads=$2           ; shift 2 ;;
    -h | --help)              usage ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

echo "library_id                : $library_id"
echo "inputbam                  : $inputbam"
echo "outputdir                 : $outputdir"
echo "outputbam                 : $outputbam"
echo "interval                  : $interval"
echo "modality                  : $modality"
echo "barcodes                  : $barcodes"
echo "validate                  : $validate"
echo "threads                   : $threads"
echo "Parameters remaining are  : $@"

# check input files
if [ ! -f $inputbam ]; then { echo "Input bam file not found"; exit 1; }; fi
if [ ! -f $barcodes ]; then { echo "Barcodes file not found"; exit 1; }; fi

# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# activate gatk conda environ
source activate gatk 

workdir=$SCRATCH1/bcfilter_${modality}/$library_id
rm -rf $workdir
mkdir -p $outputdir
mkdir -p $workdir

if [ $interval ]; then
  echo "Filtering input bam by interval: $interval"
  samtools view -bS $inputbam $interval |pv > $workdir/$interval.bam
  samtools index -@ $threads $workdir/$interval.bam
  inputbam=$workdir/$interval.bam
fi

# print out barcode and sample name "orig.ident" columns and then filter by sample name
awk -F ',' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}}{ print $(f["barcode"]),$(f["orig.ident"]) }' $barcodes|\
awk -v a="${library_id}" '{if($2 == a) {print $1"-1"}}' > /tmp/${modality}_barcodes.$library_id.txt

# filter for sample barcodes that pass QC
echo "Filtering bam file by CB barcode tag"
rm $outputdir/$outputbam 2> /dev/null
export TMPDIR=$workdir/TMPDIR  
rm -rf $TMPDIR; mkdir -p $TMPDIR 2> /dev/null
cores=$(($threads / 2))
subset-bam \
  --bam $inputbam \
  --cell-barcodes /tmp/${modality}_barcodes.$library_id.txt \
  --out-bam $outputdir/$outputbam \
  --cores $cores \
  | pv -t

# index output bam
echo "Indexing barcode filtered bam"
samtools index -@ $threads $outputdir/$outputbam

# validate the output bam with gatk and ignore warnings which result from missing NM tags
# if this step is omitted subset-bam can produce silent errors resulting in malformed bam output
if [ $validate = "true" ]; then
  gatk ValidateSamFile -I $outputdir/$outputbam --IGNORE_WARNINGS true
fi

# clean up
rm -rf $workdir

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

echo -e "\033[35;40mBarcode filtering completed in $(format_time $SECONDS)\033[0m"
