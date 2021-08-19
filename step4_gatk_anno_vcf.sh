#!/bin/bash
# this script will annotate vcf using gatk funcotator

# Set some default values:
threads=1
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

EOF
exit 1
}


if [[ ${#} -eq 0 ]]; then
   usage
fi

PARSED_ARGUMENTS=$(getopt -a -n step4_gatk_anno_vcf.sh \
-o n:v:d:o:r:sirVt:h \
--long library_id:,inputvcf:,outputdir:,outputvcf:,reference:,output_table:,modality:,funcotation:,verbose,threads:,help -- "$@")

echo "PARSED_ARGUMENTS are $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -n | --library_id)        library_id=$2        ; shift 2 ;;
    -v | --inputvcf)          inputvcf=$2          ; shift 2 ;;
    -d | --outputdir)         outputdir=$2         ; shift 2 ;;
    -o | --outputvcf)         outputvcf=$2         ; shift 2 ;;
    -r | --reference)         reference=$2         ; shift 2 ;;
    -a | --output_table)      output_table=$2      ; shift 2 ;;
    -m | --modality)          modality=$2          ; shift 2 ;;
    -f | --funcotation)       funcotation=$2       ; shift 2 ;;
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
echo "reference                 : $reference"
echo "output_table              : $output_table"
echo "modality                  : $modality"
echo "funcotation               : $funcotation"
echo "verbose                   : $verbose"
echo "threads                   : $threads"
echo "Parameters remaining are  : $@"

# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# set exit status temporary file for monitoring return values in parallel processes
echo $exit_status > /tmp/exit_status.txt

# activate gatk conda environ
source activate gatk

# load workdir
workdir=$SCRATCH1/gatk_genotype/$modality/$library_id/annotation
mkdir -p $workdir 2> /dev/null

# stream GATK output to terminal o/w capture in log file
if [ $verbose = "true" ]; then
  outputlog=/dev/stdout
elif [ $verbose = "false" ]; then
  outputlog=$workdir/log.out
fi

# create funcotation director
mkdir $outputdir 2> /dev/null

# download funcotator resource (if not already there)
if [ -d $funcotation ]; then
  echo "Funcotator resources detected"
else
  echo "Downloading funcotator resources"
  gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download --overwrite-output-file \
    -O $reference/funcotator    
fi

# update output vcf file names
funcotated_vcf=$outputdir/$outputvcf
funcotated_table=$outputdir/$output_table
mkdir $outputdir 2> /dev/null

# index vcf
gatk IndexFeatureFile -I $inputvcf

# create scatter gather intervals across no. threads
gatk SplitIntervals \
  -R $reference/fasta/genome.fa \
  -O /tmp/interval_files_folder \
  -L $inputvcf \
  --scatter-count $threads

# # create array of scatter gather intervals  
scatter_intervals=$(ls /tmp/interval_files_folder)

# remove vcf list if session has been used multiple times
if [ -f /tmp/vcf.list ]; then
  rm /tmp/vcf.list
fi

echo "Annotating contigs with funcotator"
for scatter_interval in ${scatter_intervals[@]}; do
  echo "$workdir/$library_id.$modality.$scatter_interval.funcotated.vcf.gz" >> /tmp/vcf.list
  # annotate variants
  gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=1" Funcotator \
    --variant $inputvcf \
    --reference $reference/fasta/genome.fa \
    --ref-version hg38 \
    --data-sources-path $funcotation  \
    --output $workdir/$library_id.$modality.$scatter_interval.funcotated.vcf.gz \
    --output-file-format VCF \
    -L /tmp/interval_files_folder/$scatter_interval \
    --disable-sequence-dictionary-validation true \
    --verbosity INFO >> ${outputlog} 2>&1 \
    || { echo -e "\033[0;33mFuncotator failed on $scatter_interval. Check $outputlog for additional info\033[0m"; echo 1 > /tmp/exit_status.txt; } &
done | pv -t
 # exit with 1 if interval failed
wait
exit_status=$(head -n1 /tmp/exit_status.txt)
if [ $exit_status -eq 1 ]; then exit 1; fi

# merge and index vcfs. GatherVcfs throws an unexpected error here that may be due to GATK not correctly ordering the contigs after splitintervals
# or multiple variants with the same context in different intervals
gatk MergeVcfs -I /tmp/vcf.list -O $funcotated_vcf

# index merged vcf
echo "Indexing merged vcf"
bcftools index --threads $threads $funcotated_vcf

echo "Processing funcotated vcf"
bcftools query -f'[%CHROM,%POS,%REF,%ALT,%GT,%FILTER\n]' $funcotated_vcf > $workdir/variants_no_header.csv
echo "CHROM,POS,REF,ALT,GT,FILTER" | cat - $workdir/variants_no_header.csv > $workdir/variants.csv

# this will parse a vcf annotated by funcotator and create headers
echo "Extracting FUNCOTATION headers from INFO field"
bcftools view $funcotated_vcf|grep ID=FUNCOTATION|cut -d: -f2|sed s'/ //g'|sed s'/">//g'|sed s'/]//g'|cut -d'|' -f1-130|sed 's/|/,/g' > $workdir/funco_headers.csv

# grab the funcotation switch to csv and and cut first 130 fields
echo "Converting FUNCOTATION in INFO field to csv"
bcftools query -f'[%INFO/FUNCOTATION\n]' $funcotated_vcf|cut -d',' -f1|cut -d'|' -f1-129|sed 's/|/,/g'|sed s'/\[//g' > $workdir/funco_no_header.csv

# merge to csv file
echo "Merging headers with annotation"
cat $workdir/funco_headers.csv $workdir/funco_no_header.csv > $workdir/funcotation.csv

echo "Generating csv FUNCOTATION file"
paste -d',' $workdir/variants.csv $workdir/funcotation.csv > $workdir/formatted.variant.csv

# filter for desired annotation columns and add a "." to empty fields,
echo "Filtering for specified column annotations and writing file"
awk -v cols='CHROM,POS,REF,ALT,GT,FILTER,Gencode_27_variantClassification,Gencode_27_codonChange,gnomAD_exome_AF,gnomAD_genome_AF,\
Gencode_27_transcriptExon,Gencode_27_hugoSymbol' \
'BEGIN{FS=OFS=","; nc=split(cols, a, ",")} NR==1{for (i=1; i<=NF; i++) hdr[$i]=i} \
{for (i=1; i<=nc; i++) if (a[i] in hdr) printf "%s%s", $hdr[a[i]], (i<nc?OFS:ORS)}' $workdir/formatted.variant.csv > $workdir/funco_no_varid.csv

# add a variant_id column
awk 'BEGIN{FS=OFS=","} {print (NR>1?$1"_"$2"_"$3"_"$4:"variant_id"), $0}' $workdir/funco_no_varid.csv > $funcotated_table


if [ $exit_status -eq 0 ]; then
  rm -rf $workdir
  format_time() {
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "%02d:%02d:%02d\n" $h $m $s
   }

  echo -e "\e[92mWriting $outputvcf to $outputdir\033[0m"

  echo -e "\033[35;40mAnnotation completed in $(format_time $SECONDS)\033[0m"
fi