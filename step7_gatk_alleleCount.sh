#!/bin/bash
# allele specific expression in single cell data using gatk4 docker image
# bam files are obtained from the wasp pipeline

# Set some default values:
inputvcf=""
waspbam=""
barcodes=""
genotype=""
library_id=""
modality=""
interval=""
pseudobulk_counts=FALSE
celltype_pseudobulk_counts=FALSE
sc_counts=FALSE
isphased=FALSE
threads=1

function usage {
        echo "Usage: $(basename $0) [-viognmlCcspt]" 2>&1
        echo '   -v  | --inputvcf           STR   path/to/input.vcf.gz eg. [vcfdir/funcotation/Control_1.pass.joint.hcphase.funco.vcf.gz]'
        echo '   -i  | --inputbam           STR   path/to/wasp.bam eg. [project/wasp_rna/Control_1.phase.wasp.bam]'
        echo '   -o  | --outputdir          STR   path/to/output directory eg. [project/wasp_rna/counts]'
        echo '   -b  | --barcodes           STR   path/to/barcodes.csv eg. [barcodes/rna_barcodes.csv]'
        echo '   -g  | --genotype           STR   genotype: [rna] [atac] [joint]'
        echo '   -n  | --library_id         STR   library_id: eg. [Control_1]'
        echo '   -m  | --modality           STR   sequencing modality for short variant discovery: [rna] [atac]'
        echo '   -l  | --interval           STR   optional: count a specified interval eg. [chr10]'
        echo '   -C  | --pseudobulk_counts        allele-specific counts with all cells grouped together'
        echo '   -c  | --celltype_counts          allele-specific counts after grouping cells by barcode annotation'
        echo '   -s  | --single_cell_counts       single cell allele-specific counts for provided barcodes'        
        echo '   -p  | --isphased                 optional: input vcf is phased. Default=[false]'
        echo '   -t  | --threads            INT   number of threads. Default=[1]'
        echo '   -h  | --help                     show usage'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

PARSED_ARGUMENTS=$(getopt -a -n step7_wasp.sh -o v:i:o:b:g:n:m:l:Ccspt:h --long inputvcf:,inputbam:,outputdir:,barcodes:,genotype:,library_id:,modality:,interval:,\
pseudobulk_counts,celltype_counts,single_cell_counts,isphased,threads:,help -- "$@")

echo "PARSED_ARGUMENTS are $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -v | --inputvcf)            inputvcf=$2                     ; shift 2 ;;
    -i | --inputbam)            waspbam=$2                      ; shift 2 ;;
    -o | --outputdir)           outputdir=$2                    ; shift 2 ;;
    -b | --barcodes)            barcodes=$2                     ; shift 2 ;;
    -g | --genotype)            genotype=$2                     ; shift 2 ;;
    -n | --library_id)          library_id=$2                   ; shift 2 ;;
    -m | --modality)            modality=$2                     ; shift 2 ;;
    -l | --interval)            interval=$2                     ; shift 2 ;;
    -C | --pseudobulk_counts)   pseudobulk_counts=TRUE          ; shift 1 ;;
    -c | --celltype_counts)     celltype_pseudobulk_counts=TRUE ; shift 1 ;;
    -s | --single_cell_counts)  sc_counts=TRUE                  ; shift 1 ;;
    -p | --isphased)            isphased=TRUE                   ; shift 1 ;;
    -t | --threads)             threads=$2                      ; shift 2 ;;
    -h | --help)                usage ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

echo "library_id                : $library_id"
echo "inputbam                  : $waspbam"
echo "outputdir                 : $outputdir"
echo "barcodes                  : $barcodes"
echo "genotype                  : $genotype"
echo "inputvcf                  : $inputvcf"
echo "modality                  : $modality"
echo "interval                  : $interval"
echo "pseudobulk_counts         : $pseudobulk_counts"
echo "celltype_counts           : $celltype_pseudobulk_counts"
echo "single_cell_counts        : $sc_counts"
echo "isphased                  : $isphased"
echo "threads                   : $threads"
echo "Parameters remaining are  : $@"


# ensure gatk and miniconda are in path when working in LSF environment
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# activate gatk conda environ
source activate gatk

# prepare a fasta dict file using the cellranger ref
# gatk CreateSequenceDictionary -R /ref/fasta/genome.fa
scbamdir=$outputdir/scbam/$library_id
workdir=$SCRATCH1/wasp_${modality}/$genotype/$library_id
mkdir -p $outputdir
mkdir -p $workdir

# add filter tag in format field to heterozygous variants
echo "Adding heterozygosity filter to input VCF"
gatk VariantFiltration \
-V $inputvcf \
-O /tmp/isHet.$(basename $inputvcf) \
--genotype-filter-expression "isHet == 1" \
--genotype-filter-name "isHetFilter"

echo "Filtering multiallelic and non-heterozygous variants from $inputvcf and retaining biallelic SNV"
# collapse variants with the same context and remove multiallelic variants
bcftools norm /tmp/isHet.$(basename $inputvcf) -m +snps |bcftools view -Oz -m2 -M2 -v snps > /tmp/single_context.vcf.gz

# remove variants that are not heterozygous
(bcftools view -h /tmp/single_context.vcf.gz; bcftools view -H /tmp/single_context.vcf.gz|grep 'isHetFilter')|\
  bcftools view -Oz - > /tmp/filter.$(basename $inputvcf)
gatk IndexFeatureFile -I /tmp/filter.$(basename $inputvcf)

# limit to specified intervals for testing purposes if -L flag selected
# specify intervals to eliminate alt contigs
# update output vcf file names
rm -rf /tmp/interval_files_folder 2> /dev/null
if [ $interval ]; then
  echo "Selected interval is $interval"
  workdir=$SCRATCH1/wasp_${modality}/$genotype/$library_id/$interval
  mkdir -p $workdir

  # filter cellranger bam for selected interval
  echo "Filtering bam for selected interval $interval"
  samtools view -@ $threads -bS $waspbam $interval|pv|samtools sort -@ $threads > $workdir/$(basename $waspbam .wasp.bam).${interval}wasp.bam
  waspbam=$workdir/$(basename $waspbam .wasp.bam).${interval}wasp.bam
  samtools index -@ $threads $waspbam
  # create scatter gather intervals across no. threads
  gatk SplitIntervals \
  -R ${modality}_ref/fasta/genome.fa \
  -O /tmp/interval_files_folder \
  --scatter-count $threads \
  -L $interval \
  --interval-set-rule INTERSECTION
else
  # create scatter gather intervals across no. threads
  gatk SplitIntervals \
  -R ${modality}_ref/fasta/genome.fa \
  -O /tmp/interval_files_folder \
  --scatter-count $threads
fi
# create array of split intervals
split_intervals=$(ls /tmp/interval_files_folder)

#########################################################
##################PSEUDOBULK COUNTS######################
# ase counting for all cell types grouped together
if [ $pseudobulk_counts == "TRUE" ]; then
  # create a new directory for split tables
  rm -rf /tmp/gather_tables; mkdir /tmp/gather_tables > /dev/null

  # do allele counting across intervals
  for split_interval in ${split_intervals[@]}; do
    # perform allele specific counting for all celltypes (pseudobulk)
    # tool will automatically filter out non-heterozygous positions and duplicates (if not already removed)
    gatk ASEReadCounter \
    -R ${modality}_ref/fasta/genome.fa \
    -I $waspbam \
    -V /tmp/filter.$(basename $inputvcf) \
    -L /tmp/interval_files_folder/$split_interval \
    -O /tmp/gather_tables/$split_interval &
  done
  wait

  echo "Gathering $celltype split interval count tables"
  # gather all the count tables into a single table
tables=($(ls /tmp/gather_tables))
(head -n1 /tmp/gather_tables/$tables; tail -q -n+2 /tmp/gather_tables/*) > $outputdir/$library_id.${interval}counts.pseudobulk.table

# add phased genotypes to counts table if input vcf is phased
  if [ $isphased == "TRUE" ]; then
    # create a genotype table with unique variant ids and their corresponding genotypes using the filtered input vcf
    bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' /tmp/filter.$(basename $inputvcf) |\
    awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2"_"$3"_"$4, $5}' > /tmp/genotype_table.tsv

    # create a unique identifier for each variant in the count table and put it in the first column
    count_table=$outputdir/$library_id.${interval}counts.pseudobulk.table
    awk 'BEGIN{FS=OFS="\t"} {print (NR>1?$1"_"$2"_"$4"_"$5:"variant_id"), $0}' $count_table > /tmp/variantid_count_table.tsv

    # join the count table and genotype table by unique variant_id in the first column
    join  -j 1 -t $'\t' <(sort /tmp/genotype_table.tsv) <(sort /tmp/variantid_count_table.tsv) > /tmp/phased_count.table

    # reheader the merged table
    (echo -e "variant_id\tGT\t$(head -n1 $count_table)"; cat /tmp/phased_count.table) > $outputdir/$library_id.${interval}counts.pseudobulk.phased.table
  fi

  # TODO: annotate count table with gnomad AF and gene names etc from FUNCOTATION

fi 
##################################################################
##################CELLTYPE PSEUDOBULK COUNTS######################
# ase counting for individual celltypes (ie group all cells of the same type and count together)
# create another bam file that only contains reads for a specified cell type
# read in the barcodes with the following format:barcode,celltype,lowres.celltype,orig.ident
if [ $celltype_pseudobulk_counts == "TRUE" ]; then
  echo "Performing cell type pseudobulk counts"
  # first format the barcode file by printing out barcode,orig.ident,lowres.celltype columns and then filter by sample name
  awk -F ',' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}}{ print $(f["barcode"]),$(f["orig.ident"]),$(f["celltype"]) }' $barcodes |\
  awk -v a="${library_id}" '{if($2 == a) {print "CB:Z:"$1"-1",$2,$3}}' > /tmp/${modality}_barcodes.$library_id   

  celltype_groups=$(tail -n+2 /tmp/${modality}_barcodes.$library_id|awk '{print $3}'|sort|uniq)
  for celltype in $celltype_groups; do
    # filter barcodes for selected celltype barcodes
    echo $celltype
    awk -v a="${celltype}" '$3 == a' /tmp/${modality}_barcodes.$library_id | \
    awk '{print $1}' | cut -d ':' -f3 > /tmp/${modality}_barcodes.$library_id.$celltype.txt 

    rm $workdir/$library_id.$celltype.${interval}wasp.bam 2> /dev/null
    TMPDIR=$workdir/TMPDIR  
    rm -rf $TMPDIR; mkdir -p $TMPDIR 2> /dev/null
    subset-bam \
    --bam $waspbam \
    --cell-barcodes /tmp/${modality}_barcodes.$library_id.$celltype.txt  \
    --out-bam $workdir/$library_id.$celltype.${interval}wasp.bam \
    --cores $threads

    # index the celltype bam to enable traversal by intervals
    samtools index -@ $threads $workdir/$library_id.$celltype.${interval}wasp.bam

    # create a new directory for split tables
    rm -rf /tmp/gather_tables; mkdir /tmp/gather_tables > /dev/null

    # do allele counting across intervals
    for split_interval in ${split_intervals[@]}; do
      # perform allele specific counting for all celltypes (pseudobulk)
      # tool will automatically filter out non-heterozygous positions and duplicates (if not already removed)
      gatk ASEReadCounter \
      -R ${modality}_ref/fasta/genome.fa \
      -I $workdir/$library_id.$celltype.${interval}wasp.bam \
      -V /tmp/filter.$(basename $inputvcf) \
      -L /tmp/interval_files_folder/$split_interval \
      -O /tmp/gather_tables/$split_interval &
    done
    wait

    echo "Gathering celltype pseudobulk count tables"
    # gather all the count tables into a single table
    tables=($(ls /tmp/gather_tables))
    (head -n1 /tmp/gather_tables/$tables; tail -q -n+2 /tmp/gather_tables/*) > $outputdir/$library_id.${interval}counts.$celltype.table

    # add phased genotypes to counts table if input vcf is phased
    if [ $isphased == "TRUE" ]; then
      # create a genotype table with unique variant ids and their corresponding genotypes using the filtered input vcf
      bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' /tmp/filter.$(basename $inputvcf) |\
      awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2"_"$3"_"$4, $5}' > /tmp/genotype_table.tsv

      # create a unique identifier for each variant in the count table and put it in the first column
      count_table=$outputdir/$library_id.${interval}counts.$celltype.table
      awk 'BEGIN{FS=OFS="\t"} {print (NR>1?$1"_"$2"_"$4"_"$5:"variant_id"), $0}' $count_table > /tmp/variantid_count_table.tsv

      # join the count table and genotype table by unique variant_id in the first column
      join  -j 1 -t $'\t' <(sort /tmp/genotype_table.tsv) <(sort /tmp/variantid_count_table.tsv) > /tmp/phased_count.table

      # reheader the merged table
      (echo -e "variant_id\tGT\t$(head -n1 $count_table)"; cat /tmp/phased_count.table) > $outputdir/$library_id.${interval}counts.$celltype.phased.table
    fi

    # TODO: annotate count table with gnomad AF and gene names etc from FUNCOTATION

  done | pv -t > /dev/null
fi
###################################################################
##################SINGLE CELL BAMS AND COUNTS######################
# generate individual single cell bams and perform single cell counts
if [ $sc_counts == "TRUE" ]; then
  echo "Creating single cell bams and performing single cell counts"

  # create output subdirectories for single cell bams
  rm -rf $scbamdir; mkdir -p $scbamdir/{counts,bam}

  # # convert vcf to bed
  bcftools query -f '%CHROM\t%POS\t%POS\n' /tmp/filter.$(basename $inputvcf) > /tmp/$(basename $inputvcf .vcf.gz).bed

  # # filter by region
  echo "Filtering wasp bam by vcf sites"
  samtools view -bS -L /tmp/$(basename $inputvcf .vcf.gz).bed $waspbam | pv > $scbamdir/$library_id.${interval}sites.bam
  bamsites=$scbamdir/$library_id.${interval}sites.bam
  samtools index -@ $threads $bamsites

  # make list of all unique barcodes in wasp bam
  echo "Retrieving unique barcodes from $(echo $(basename $bamsites))"
  barcodes=($(samtools view $bamsites | pv | cut -f 12- | tr "\t" "\n"  | grep  "^CB:Z:"  | cut -d ':' -f3 | sort | uniq )) 
  num_barcodes=${#barcodes[@]}

  # split wasp bam file into cell-specific bams in parallel loop using subset-bam
  echo "Splitting $(echo $(basename $bamsites)) into single cell bam files"
  TMPDIR=$workdir/TMPDIR
  rm -rf $TMPDIR; mkdir -p $TMPDIR 2> /dev/null
  for barcode in ${barcodes[*]}; do \
    echo $barcode > /tmp/barcode.txt
    subset-bam \
    --bam $bamsites \
    --cell-barcodes /tmp/barcode.txt \
    --out-bam $scbamdir/bam/$barcode.bam \
    --cores $threads
    echo $barcode 
  done | pv -l -s $num_barcodes > /dev/null

  # allele specific counts of celltype sam files with gatk in parallel loop
  # output as [barcode].counts in counts dir
  echo "Counting single cell bam files"
  bamfiles=($(ls $scbamdir/bam))
  for bamfile in ${bamfiles[*]}; do
    table=$(basename $bamfile .bam).counts
    gatk ASEReadCounter \
      -R ${modality}_ref/fasta/genome.fa \
      -I $scbamdir/bam/$bamfile \
      -V /tmp/filter.$(basename $inputvcf) \
      --verbosity ERROR \
      -O $scbamdir/counts/$table > /dev/null 2>&1  &
    while (( $(jobs |wc -l) >= (( ${threads} + 1 )) )); do
  	sleep 0.1
    done
    echo $bamfile
  done | pv -l -s $num_barcodes > /dev/null	 

  # add a barcode column to each count table
  tables=($(ls $scbamdir/counts))
  rm -rf /tmp/gather_tables; mkdir /tmp/gather_tables > /dev/null
  for table in ${tables[@]}; do
    barcode=$(basename $table .counts)
    awk -v d="$barcode" 'BEGIN{FS=OFS="\t"} {print $0, (NR>1?d:"barcode")}' $scbamdir/counts/$table > /tmp/gather_tables/$table &
  done 

  # concatenate all the barcode counts into a single table
  (head -n1 /tmp/gather_tables/$tables; tail -q -n+2 /tmp/gather_tables/*) > $outputdir/$library_id.${interval}counts.single_cell.table

  # add phased genotypes to counts table if input vcf is phased
  if [ $isphased == "TRUE" ]; then
    # create a genotype table with unique variant ids and their corresponding genotypes using the filtered input vcf
    bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' /tmp/filter.$(basename $inputvcf) |\
    awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2"_"$3"_"$4, $5}' > /tmp/genotype_table.tsv

    # create a unique identifier for each variant in the count table and put it in the first column
    count_table=$outputdir/$library_id.${interval}counts.single_cell.table
    awk 'BEGIN{FS=OFS="\t"} {print (NR>1?$1"_"$2"_"$4"_"$5:"variant_id"), $0}' $count_table > /tmp/variantid_count_table.tsv

    # join the count table and genotype table by unique variant_id in the first column
    join  -j 1 -t $'\t' <(sort /tmp/genotype_table.tsv) <(sort /tmp/variantid_count_table.tsv) > /tmp/phased_count.table

    # reheader the merged table
    (echo -e "variant_id\tGT\t$(head -n1 $count_table)"; cat /tmp/phased_count.table) > $outputdir/$library_id.${interval}counts.single_cell.phased.table
  fi 

  # TODO: annotate count table with gnomad AF and gene names etc from FUNCOTATION

  # clean up single cell bam directory and individual count tables
  rm -rf $scbamdir
fi
#########################################################
#########################################################

# cleanup
rm -rf $workdir

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

echo -e "\033[35;40mAllele counting completed in $(format_time $SECONDS)\033[0m"