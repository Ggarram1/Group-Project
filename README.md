# Diversity Among Homologous Genes #  

## Description ## 

Many papers use mice and other organisms as comparative models for human-centric questions. A paper published in Nature uses comparative transcriptomics to assess the homology across species at various levels. In figure 2 of the paper, transcriptional profiles are used to cluster samples based on individual homology (same species, different organs) and species homology (different species same organ). For our project, we will recreate this comparison using the available data. If possible, we may attempt to use available transcriptomes to analyze an additional species.

## Example Published Figure ## 
[Example Published Figure](https://www.nature.com/articles/s41586-019-1338-5/figures/1) Our project will focus on the data presented in panel B. 

## Data set ## 
[Datasets with links](https://apps.kaessmannlab.org/evodevoapp/) We will retrieve the datasets from this website to track the differences in the transcriptomes between species. 
[Datasets with links](https://useast.ensembl.org/index.html) We will use ensemble to identify homologs across different species in order to utilize the data effectively. 

## Software ## 
~~[Software with links](https://www.python.org) We will utilize python when deciphering and ordering our datasets to each other.~~
- We have found a more easily digestible way to organize our datasets still in R! This will help us in maintaing the same language, as well as keeping our data within the same software. 
- [Software with links](https://www.r-project.org) We will utilize R script when plotting the principle components between species transcriptomes.

## Proposed Steps ## 

1) Download data for organisms ✅
2) ~~Plot the data to recreate the figure above on a smaller scale - for a single organism - to familiarize ourselves with PCA~~ 
2) Establish method using ensembl to find homologs  ✅
3) Join datasets together ✅
4) Organize data - sort out how to handle both tissue/age/replicate as well as homolog replicates ✅
 --> folllowing recommendations, we are using the tissue/age/replicate details as metadata and only using 1:1 homologs
5) First use two datasets together for proof of concept, perform VST and plot PCA ✅
6) Append all different organismal data to make one master file ⏳
7) Plot the full data!


## Presentation Link ##
https://docs.google.com/presentation/d/134gtdcx9AfGN8t4SqgLmlcO4cLPa38ICfstJ2mVbg4k/edit#slide=id.g319c96af645_0_46
