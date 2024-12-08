


#load libraries
library(biomaRt)
library(tidyverse)
library(dplyr)

#load in human gene data from Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#define function to pull out homologs
get_homologs <- function(gene_ids, source_dataset, target_dataset) {
  #connect to Ensembl Biomart
  ensembl <- useEnsembl(biomart = "genes", dataset = source_dataset)
  
  #retrieve homologs (human to mouse)
  homologs <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 
                   'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name'),
    filters = 'ensembl_gene_id',
    values = gene_ids,
    mart = ensembl
  )
  
  return(homologs)
}

#load in human mouse expression data
Human_Exp <- read_delim("Human_rpkm.txt")
Mouse_Exp <- read_delim("Mouse_rpkm.txt")

#get human gene IDs from Human_Exp
human_gene_ids <- Human_Exp$Names 

#get homologs for human genes and mouse genes
human_homologs <- get_homologs(human_gene_ids, 'hsapiens_gene_ensembl', 'mmusculus_gene_ensembl')
head(human_homologs)

#rename for clarity
Mouse_Exp <- Mouse_Exp %>%
  rename(ensembl_mm_gene_id = Names)

#left join Mouse_Exp with human homologs to get corresponding human Ensembl gene IDs
Mouse_human <- left_join(Mouse_Exp, human_homologs, 
                         by = c("ensembl_mm_gene_id" = "mmusculus_homolog_ensembl_gene"))

#mutate the column to use human gene IDs where available
Mouse_human <- Mouse_human %>%
  mutate(ensembl_mm_gene_id = ifelse(!is.na(ensembl_gene_id), ensembl_gene_id, ensembl_mm_gene_id))

#check merged data
head(Mouse_human)



