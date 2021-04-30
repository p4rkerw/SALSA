# **Welcome to 🌶️ SALSA! - Single Cell Allele-Specific Analysis**

Here you will find a workflow to analyze single cell allele-specific expression and chromatin accessibility obtained from 10X Genomics Single Cell Gene Expression and Single Cell ATAC datasets. To get started follow along with the [tutorial](https://github.com/p4rkerw/SALSA/tree/main/Tutorial/README.md). 

If you use 🌶️SALSA please cite the following manuscript:
```
Muto Y, Wilson PC, Ledru N, Wu H, Dimke H, Waikar SS, Humphreys BD.<br/> Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney. Nat Commun. 2021 Apr 13;12(1):2190. doi: 10.1038/s41467-021-22368-w. PMID: 33850129.
```
Thanks!  
Parker
<br/><br/>
![alt text](http://humphreyslab.com/wp-content/uploads/2015/12/favicon-H.jpg)  
Visit the Humphrey's lab website:   
www.humphreyslab.com  
<br/>
Check out our interactive datasets with Kidney Interactive mulTiomics (KIT):  
http://humphreyslab.com/SingleCell/
<br/><br/>
Find us on Twitter: 
<br/>
  <a href="https://twitter.com/parkercwilson?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> @parkercwilson</a>
  <a href="https://twitter.com/HumphreysLab?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> @HumphreysLab</a>
<br/><br/>
Find us on Docker Hub:  
[p4rkerw@dockerhub](https://hub.docker.com/search?q=p4rkerw&type=image)
<br/>

**Allele Specific Analysis:**    
These scripts can be run in a publicly-available docker container found here: [p4rkerw@dockerhub](https://hub.docker.com/search?q=p4rkerw&type=image)

(Follow the steps in order) 
1. Genotype snRNA or snATAC libraries using GATK
2. (Optional) - Merge genotypes obtained from matched snRNA-snATAC or multimodal libraries
3. (Recommended) - Phase the genotype using shapeit and the 1000G reference    
4. (Optional) - Annotate the vcf with GATK Funcotator to evaluate gnomAD MAF and variant context  
5. Filter the 10X genomics positions sorted bam file by cell barcode using the subsetbam utility  
6. Apply the WASP pipeline to barcode-filtered bam files to perform variant-aware realignment  
7. Obtain pseudobulk, cell-specific, or single cell allele-specific counts with GATK ASEReadCounter  

