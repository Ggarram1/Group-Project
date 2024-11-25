# Diversity Among Homologous Genes #  
## Prior Feedback ##

Would it be possible to limit your analysis to only the 1-to-1 homologs, which are more likely to be orthologs? That might lower the number of genes under consideration, but I am guessing that you would still have a lot! You would basically be throwing out any genes that have 1-to-many or many-to-1 patterns of homology.
- Yes this is possible! And based upon your feedback and recommendations we have narrowed our analysis to only the 1-to-1 homologs. By doing so we have disregarded any genes that had different patterns of homology and significantly reduced the number of genes we are considering for this analysis. In the current code we have set up a function which is filtering throguh each gene and counting how many times it appears then selecting for only the ones that appear once. 


I think this can all be essentially considered metadata. So the entire gene expression matrix gets input to PCA, but then you can do things like color your datapoints or set the shapes of your datapoints in your figure based on these various metadata features and answer questions about whether species versus tissue versus replicate versus age is more important for explaining gene expression similarity/differences.
- We have not completed the code for the PCA/dimension reduction but moving forward we will be taking this metadata approach when we plot the values using ggplot. Now that we have the correct 1-to-1 homologs the analysis should run smoothly from this point. When ploting our data we can split up the various features and assess which feature has a more predominant distribution, thus explanation some of the variation we observe. Depending on the difficulty of the analysis from this point we can deliberate whether we can apply these methods to compare other species and their homologs. 



## New progress since last submission ##
- Since the prior submission we have written a code (and uploaded it as Filtered Homologs) that allows us to filter our viable gene data by removing any gene homologs that are not 1-to-1. This code can be seen in lines 52-72 where we select our unique human gene data and filtering through each of the genes (columns). It will summarize each time the gene appears and then filter only for the genes that appear once thus disregarding any 1-to-many or vice-versa. It will then do the same for the mouse genes and the final lines of code will look at the intersection of unique mouse and human homologs to find those which appear in both organisms only one time. 
- We have also now performed the VST and an initial PCA on the mouse and human data. This first PCA appears quite different to the published one, but that is to be expected. We are starting with two higher mammals so, the PCA will return different results than when analyzing a greater range of organisms

## Struggles you are encountering and questions you would like advice on ##
~~We had many struggles in joining the two datasets but once we got advice from Sadhana and Fred we were able to fix this issue and ultilize the left join function. However at this point, we don't have any major questions and we haven't had enough time to find a new struggle.~~
We have made quite a bit of progress on the trial two organism PCA. Given the current code we uploaded, is there a less messy way to organize everything? It is getting a bit hard to follow. Also, the age key is very ugly. How can we make it look more like the one in the paper?
