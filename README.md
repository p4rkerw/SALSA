# SALSA
Single Cell Allele-Specific Analysis

Welcome to our GitHub repository!  
Here you will find a workflow and models for analyzing single cell allele-specific expression obtained from 10X Genomics scRNA-seq datasets.  
<br/>
Thanks,  
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
These scripts can be run in publicly-available docker containers found here: [p4rkerw@dockerhub](https://hub.docker.com/search?q=p4rkerw&type=image)

Each script has an example command to run the corresponding docker container  

(Follow the steps in order) 
1. Genotype the snRNA or snATAC libraries using GATK (or obtain a vcf from another method)
2. (Optional) - Merge genotypes obtained from matched snRNA-snATAC or multimodal libraries
3. Phase the genotype using shapeit and the 1000G reference    
4. (Optional) - Annotate the vcf with GATK Funcotator to evaluate gnomAD MAF and variant context  
5. Filter the 10X genomics positions sorted bam file by cell barcode using the subsetbam utility  
6. Apply the WASP pipeline to barcode-filtered bam files to perform variant-aware realignment  
9. Obtain pseudobulk, cell-specific, or single cell allele-specific counts with GATK ASEReadCounter  
10. Analyze allele-specific expression across multiple samples with ASEP  

