#load libraries
library(biomaRt)
library(tidyverse)
library(dplyr)

#load in human gene data from Ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#define function to pull out homologs

#The attributes need to change, theres a parameter that isnt being used
# Go through biomart and find the other attributes for the other organisms
get_homologs <- function(gene_ids, source_dataset, target_one, target_two) {
  #connect to Ensembl Biomart
  ensembl <- useEnsembl(biomart = "genes", dataset = source_dataset)
  
  #retrieve homologs (human to mouse)
  homologs <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', target_one, target_two), 
    filters = 'ensembl_gene_id',
    values = gene_ids,
    mart = ensembl
  )
  
  return(homologs)
}

#load in human mouse expression data
Human_Exp <- read_delim("Human_rpkm.txt")
Mouse_Exp <- read_delim("Mouse_rpkm.txt")
Chicken_Exp <- read_delim("Chicken_rpkm.txt")
Macaque_Exp <- read_delim("Macaque_rpkm.txt")
Opossum_Exp <- read_delim("Opossum_rpkm.txt")

#get human gene IDs from Human_Exp
human_gene_ids <- Human_Exp$Names 

#get homologs for human genes and other species genes
human_to_mouse_homologs <- get_homologs(human_gene_ids, 'hsapiens_gene_ensembl',
                                        'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name')
human_to_macaque_homologs <- get_homologs(human_gene_ids, 'hsapiens_gene_ensembl', 
                                          'mmulatta_homolog_ensembl_gene', 'mmulatta_homolog_associated_gene_name')
human_to_opossum_homologs <- get_homologs(human_gene_ids, 'hsapiens_gene_ensembl', 
                                          'mdomestica_homolog_ensembl_gene', 'mdomestica_homolog_associated_gene_name')
human_to_chicken_homologs <- get_homologs(human_gene_ids, 'hsapiens_gene_ensembl',
                                          'ggallus_homolog_ensembl_gene', 'ggallus_homolog_associated_gene_name')

#Honestly when we get to this step, lets just copy paste and change names 

#first, take the human homologs from before and pull out those which only appear
#one time using the summarize and filter functions
#now, take the mouse homologs and again only take those which appear once
Unique_human <- human_to_mouse_homologs %>% 
  dplyr::select(ensembl_gene_id) %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(n()) %>% 
  filter(`n()` == 1) %>% 
  pull(ensembl_gene_id)

Unique_mouse <- human_to_mouse_homologs %>% 
  dplyr::select(mmusculus_homolog_ensembl_gene) %>% 
  group_by(mmusculus_homolog_ensembl_gene) %>% 
  summarize(n()) %>% 
  filter(`n()` == 1) %>% 
  pull(mmusculus_homolog_ensembl_gene)

Unique_chicken <- human_to_chicken_homologs %>% 
  dplyr::select(ggallus_homolog_ensembl_gene) %>% 
  group_by(ggallus_homolog_ensembl_gene) %>% 
  summarize(n()) %>% 
  filter(`n()` == 1) %>% 
  pull(ggallus_homolog_ensembl_gene)

Unique_macaque <- human_to_macaque_homologs %>% 
  dplyr::select(mmulatta_homolog_ensembl_gene) %>% 
  group_by(mmulatta_homolog_ensembl_gene) %>% 
  summarize(n()) %>% 
  filter(`n()` == 1) %>% 
  pull(mmulatta_homolog_ensembl_gene)

Unique_opossum <- human_to_opossum_homologs %>% 
  dplyr::select(mdomestica_homolog_ensembl_gene) %>% 
  group_by(mdomestica_homolog_ensembl_gene) %>% 
  summarize(n()) %>% 
  filter(`n()` == 1) %>% 
  pull(mdomestica_homolog_ensembl_gene)

#now, look at the intersection of unique mouse and human homologs to find
#those which appear in both organisms only one time


#Add in all of the unique homologs from monkey, opossum, etc. to the unqiue_homologs 


unique_mouse_homologs <- human_to_mouse_homologs %>% 
  filter(ensembl_gene_id %in% Unique_human & mmusculus_homolog_ensembl_gene %in% Unique_mouse)

unique_chicken_homologs <- human_to_chicken_homologs %>% 
  filter(ensembl_gene_id %in% Unique_human & ggallus_homolog_ensembl_gene %in% Unique_chicken)

unique_macaque_homologs <- human_to_macaque_homologs %>% 
  filter(ensembl_gene_id %in% Unique_human & mmulatta_homolog_ensembl_gene %in% Unique_macaque)

unique_opossum_homologs <- human_to_opossum_homologs %>% 
  filter(ensembl_gene_id %in% Unique_human & mdomestica_homolog_ensembl_gene %in% Unique_opossum)

human_to_mouse_homologs <- human_to_mouse_homologs %>% distinct(ensembl_gene_id, mmusculus_homolog_ensembl_gene, .keep_all = TRUE)
human_to_chicken_homologs <- human_to_chicken_homologs %>% distinct(ensembl_gene_id, ggallus_homolog_ensembl_gene, .keep_all = TRUE)
human_to_macaque_homologs <- human_to_macaque_homologs %>% distinct(ensembl_gene_id, mmulatta_homolog_ensembl_gene, .keep_all = TRUE)
human_to_opossum_homologs <- human_to_opossum_homologs %>% distinct(ensembl_gene_id, mdomestica_homolog_ensembl_gene, .keep_all = TRUE)
# Example: Left joining multiple datasets based on human homologs (ensembl_gene_id)
datasets <- list(unique_mouse_homologs, unique_chicken_homologs, unique_macaque_homologs, unique_opossum_homologs)
# Perform the left join using Reduce()
combined_data <- Reduce(function(x, y) {
  left_join(x, y, by = "ensembl_gene_id")  # Join based on human homolog gene ID
}, datasets)

combined_data <- combined_data %>%
  filter(!is.na(mmusculus_homolog_ensembl_gene) & 
           !is.na(ggallus_homolog_ensembl_gene) & 
           !is.na(mmulatta_homolog_ensembl_gene) & 
           !is.na(mdomestica_homolog_ensembl_gene))

combined_data$ggallus_homolog_ensembl_gene <- trimws(combined_data$ggallus_homolog_ensembl_gene)
Chicken_Exp$Names <- trimws(Chicken_Exp$Names)
Chicken_combo <- combined_data %>% 
  left_join(Chicken_Exp, by = c("ggallus_homolog_ensembl_gene" = "Names"))

Combined_Exp <- combined_data %>% 
  left_join(Mouse_Exp, by = c("mmusculus_homolog_ensembl_gene" = "Names")) %>%
  left_join(Chicken_Exp, by = c("ggallus_homolog_ensembl_gene" = "Names")) %>%
  left_join(Macaque_Exp, by = c("mmulatta_homolog_ensembl_gene" = "Names")) %>%
  left_join(Opossum_Exp, by = c("mdomestica_homolog_ensembl_gene" = "Names"))
            
True_expression <- left_join(Combined_Exp, Human_Exp,
                             by = c("ensembl_gene_id" = "Names")) 
                        
#Below are my failed attempts at combining our data using left join, I finally got it to work using the code youve seen previously. 

#left join Mouse_Exp with human homologs to get corresponding human Ensembl gene IDs
#Mouse_Expression_data <- left_join(unique_mouse_homologs, Mouse_Exp,
                       #  by = c("mmusculus_homolog_ensembl_gene" = "Names"))
#Chicken_Expression_data <- left_join(unique_chicken_homologs, Chicken_Exp, 
                            # by = c("ggallus_homolog_ensembl_gene" = "Names"))
#Macaque_Expression_data <- left_join(unique_macaque_homologs, Macaque_Exp, 
                            # by = c("mmulatta_homolog_ensembl_gene" = "Names"))
#Opossum_Expression_data <- left_join(unique_opossum_homologs, Opossum_Exp, 
                                   #  by = c("mdomestica_homolog_ensembl_gene" = "Names"))
#left join Human_Exp with other homologs
#Expression_data <- left_join(Mouse_Expression_data, Human_Exp,
                       #  by = c("ensembl_gene_id" = "Names"))
#C_Expression_data <- left_join(Chicken_Expression_data, Human_Exp,
                      #       by = c("ensembl_gene_id" = "Names"))
#Ma_Expression_data <- left_join(Macaque_Expression_data, Human_Exp,
                        #       by = c("ensembl_gene_id" = "Names"))
#Op_Expression_data <- left_join(Opossum_Expression_data, Human_Exp,
                        #        by = c("ensembl_gene_id" = "Names"))


#True_Exp <- left_join(Expression_data, Ma_Expression_data,
                      #by = c("ensembl_gene_id" = "ensembl_gene_id"))


#mutate the column to use human gene IDs where available
#add the mutate function for each animal but not the human clumn cause thats what we want. 
Human_to_all <- True_expression %>%
  mutate(ensembl_mm_gene_id = ifelse(!is.na(ensembl_gene_id), ensembl_gene_id, ensembl_mm_gene_id),
         ensembl_gg_gene_id = ifelse(!is.na(ensembl_gene_id), ensembl_gene_id, ensembl_gg_gene_id),
         ensembl_md_gene_id = ifelse(!is.na(ensembl_gene_id), ensembl_gene_id, ensembl_gg_gene_id),
  )

All_organism_data <- Human_to_all %>% 
  dplyr::select(-c("external_gene_name.x", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name", 
                   "external_gene_name.y", "ggallus_homolog_ensembl_gene", "ggallus_homolog_associated_gene_name",
                   "external_gene_name.x.x", "mmulatta_homolog_ensembl_gene", "mmulatta_homolog_associated_gene_name", 
                   "external_gene_name.y.y", "mdomestica_homolog_ensembl_gene", "mdomestica_homolog_associated_gene_name"))

write_delim(All_organism_data, "FiveOrganism_expressionData.txt", delim = "\t")
