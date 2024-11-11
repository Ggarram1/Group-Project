library(biomaRt)
library(tidyverse)
library(dplyr)

#load in human gene data from Ensembl
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

#define function to pull out homologs
get_homologs <- function(gene_ids, source_dataset, target_dataset) {
  ensembl <- useEnsembl(biomart = "genes")
  ensembl <- useDataset(dataset = source_dataset, mart = ensembl)
  homologs <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 
                                   'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name'),
                    filters = 'ensembl_gene_id',
                    values = gene_ids,
                    mart = ensembl)
  return(homologs)
}
#load in human data
Human_Exp <- read_delim("Human_rpkm.txt")
#load in mouse data
Mouse_Exp <- read_delim("Mouse_rpkm.txt")

#get homologs for human genes and mouse genes
human_gene_ids <- Human_Exp$Names 
human_homologs <- get_homologs(human_gene_ids, 'hsapiens_gene_ensembl', 'mmusculus_gene_ensembl')
head(human_homologs)

#rename the column in Mouse_Exp for clarity
Mouse_Exp <- Mouse_Exp %>%
  rename(ensembl_mm_gene_id = Names)

#rename column in mouse data to be more informative
Mouse_Exp <- Mouse_Exp %>%
  rename(ensembl_mm_gene_id = Names)

#left join to add the mouse ortholog gene IDs to Mouse_Exp
# Left join Mouse_Exp with human_homologs to get the corresponding human Ensembl gene IDs
Mouse_human <- left_join(Mouse_Exp, human_homologs, 
                         by = c("ensembl_mm_gene_id" = "mmusculus_homolog_ensembl_gene"))


Mouse_human <- Mouse_human %>%
  mutate(ensembl_mm_gene_id = ifelse(!is.na(ensembl_gene_id), ensembl_gene_id, ensembl_mm_gene_id))

#check merging!
head(Mouse_human)







