ONT Methylation Read Clustering

Motivation: 
CARD brain samples are compromised of several cell types. Methylation atlas of pure cell types paper prove that there are regions of specific hypo/hyper methylation per cell type ( https://www.nature.com/articles/s41586-022-05580-6 ). 
Proof of Principle paper: 
https://www.nature.com/articles/s41588-022-01188-8


Aims
This work aims to identify cell types, or at least the number of potential cell types, within a sample based off the methylation profile of each read. By identifying reads that come from particular cell types can allow us to have a more fine grained analysis of differential methylation within each cell type in response to control/pathogenic brain samples. 
Note: diseased brains may complicate this analysis as their methylation profiles may be altered by disease. 

Methods: <br>
Parse the bam and get all the methyl sites in terms of reference cooridnates. <br>
Build an overlap graph of the CpG sites <br>
	Calculate overlap with a model that accounts for: <br>
* Measurement error 
* Information at each site ( degree of mixed methylation at each stie) <br>
* Remove isolated vertices <br>
* Keep significant overlaps constructs local ‘methylotpyes’ <br>
Use this graph to infer cell types within the sample 


