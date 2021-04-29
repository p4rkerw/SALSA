#! /bin/bash
# this script will generate a premrna reference for counting snRNAseq libraries with hg38

# download and unpack the ref
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz && \
	tar -xvzf refdata-gex-GRCh38-2020-A.tar.gz

# create custom reference that will count intronic reads
awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ print; $3="exon"; $9 = \
	gensub("(transcript_id\\s\"{0,1})([^;\"]+)(\"{0,1});", "\\1\\2_premrna\\3;", "g", $9); print; next}{print}' \
      refdata-gex-GRCh38-2020-A/genes/genes.gtf > GRCh38-2020-A.premrna.gtf

cellranger mkref \
--genome="GRCh38-2020-A.premrna" \
--fasta="refdata-gex-GRCh38-2020-A/fasta/genome.fa" \
--genes="GRCh38-2020-A.premrna.gtf" 
