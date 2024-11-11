# Diversity Among Homologous Genes #  
## Prior Feedback ##
"I think this goal can be achieved by building off of the live-coding and homework PCA exercises from the differential expression week of quant bio lab. A couple of questions:

What transformations of the gene expression count data are necessary before PCA?
    Although not a transformation, first the data will need to be filtered. As we are combining data from different organisms, it may introduce some noise/make the data harder to interprate if we include genes with no corresponding homologs. To do this, we will first filter out the genes with no data (not expression levels equal to zero, but those with no homologs in a given organism) across all organisms. Additionally, for simplicity, we will only maintain one gene in the case of numerous homologs for a given gene. Once the data have been adequetly filtered, we can proceed with the transformations. For the PCA, we will apply a variance stabilizing transformation (vst).


How will you use Ensembl to identify one-to-one orthologs? See this link as a guide: https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
    The code we have established so far has also been uploaded, but there are some details which still remain. Ensembl and BioMart allow us to find gene names of corresponding homologs between the species. Although we had initially looked into orthologs, the diversity between them makes it a bit outside the scope of this project, and we are electing to use the homologs instead. In cases of numerous homologs for a particular gene, we are still considering how to select which to keep (question below).
    
What functions will you use for "appending" the individual-species dataframes? What shared field will you use to merge each pair?"
    We will use left-join to add the homolog names, doing this to both dataframes, and then left-join again to append the individual-species data-frames. We will use this left-joined homolog name field to merge the pairs.

## New progress since last submission ##
    Since the prior submission, we have written a code that allows use to use BioMart/ensembl to create a function to find homologs between species. We then applied this to the human and mice data used in the paper, as a proof of concept. This has worked well for these and we will continue to apply this to the other organisms. We now need to append the human data to the mouse data (organize the tissue/age data - question below) and run the VST and a test PCA. 

## Struggles you are encountering and questions you would like advice on ##
- After speaking with Fred, we decided to only keep one gene in cases of numerous homologs (as mentioned). All the methods I can think of to do this will introduce some bias. Would it be best to keep the gene with the highest expression (hoping to reduce background) or would it be best to have a code to randomly pick one of the homologs? If it is the latter, how should we go about doing this?
- How would you recommend organizing the tissue/age/replicate data? Each organism has gene expression data from different organs taken at different stages of development and multiple replicates - the columns are titled as such. Should we compress these data in any way/should we group the data?