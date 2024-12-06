library(tidyverse)
library(broom)
library(stringr)
library(DESeq2)

#set working directory to where data and output are stored
setwd("~/Downloads/")

#load gene expression data
FourOrganism_expression <- read_delim("FourOrganism_expressionData.txt", 
                                      col_types = cols(
                                        ensembl_gene_id = col_character(),  #keep gene IDs as characters
                                        .default = col_double()             #default all other columns as numeric
                                      ))
#FourOrganism_expression <- FourOrganism_expression %>% select(-ensembl_mm_gene_id, -ensembl_gg_gene_id, -ensembl_md_gene_id)
#check on first few rows of the dataset
head(FourOrganism_expression)

#pull out gene IDs
gene_ids <- FourOrganism_expression$ensembl_gene_id

#pull out gene count data
count_matrix <- FourOrganism_expression %>%
  dplyr::select(-ensembl_gene_id) %>%
  as.data.frame()
summary(count_matrix)

str(count_matrix)
sapply(count_matrix[1:5], class)
min(count_matrix, na.rm = TRUE)
#0
max(count_matrix, na.rm = TRUE)
#196855.4 - not all the values are actually 0
sum(is.na(count_matrix))
#1119962 - this might be why
#first convert these values to 0
count_matrix[is.na(count_matrix)] <- 0
sum(is.na(count_matrix))
#0
sum(count_matrix == 0) 
#2287911 - this is many but since we are trying to avoid 1:many homologs, 
#we are going to remove these
sum(count_matrix != 0) 
#13157009 - we still have lots to work with!


#convert all columns to numeric using `mutate` to avoid creating a list
count_matrix_clean <- count_matrix %>%
  mutate(across(everything(), ~ as.numeric(.)))
#explicitly convert to a matrix to ensure it's numeric
count_matrix_clean <- as.matrix(count_matrix_clean)

#check the class and mode
class(count_matrix_clean)  
#matrix, array
mode(count_matrix_clean)   
#numeric

#remove rows containing 0 values
count_matrix_clean <- count_matrix %>%
  mutate(across(everything(), ~ as.numeric(.))) %>%
  as.matrix() %>%
  .[rowSums(. != 0) > 0, ]
#convert count data to integers by rounding
count_matrix_clean <- round(count_matrix_clean)

dim(count_matrix_clean) 
sum(rowSums(count_matrix_clean != 0) == 0) 

#use gene_ids as row names
rownames(count_matrix_clean) <- gene_ids[rownames(count_matrix_clean) %in% gene_ids]
count_matrix_clean <- count_matrix_clean[, !colnames(count_matrix_clean) %in% "ensembl_mm_gene_id"]
#pull out organ and age metadata from column names
organ <- str_extract(colnames(count_matrix_clean), "^[^\\.]+") 
sample_names <- colnames(count_matrix_clean)
age <- str_extract(sample_names, "(?<=\\.)[^.]+")
organism <- organism[1:ncol(count_matrix_clean)]
organism <- rep(NA, ncol(count_matrix_clean))
organism[1:316] <- "Mouse"     
organism[317:481] <- "Macaque"   
organism[482:713] <- "Possum"  
organism[714:1010] <- "Human"
organism <- factor(organism)
#organism <- case_when(
 # str_starts(sample_names, "mouse") ~ "Mouse",           # For Mouse
  #str_starts(sample_names, "human") ~ "Human",           # For Human
  #str_starts(sample_names, "opossum") ~ "Opossum",       # For Opossum
  #str_starts(sample_names, "macaque") ~ "Macaque",       # For Macaque
  #TRUE ~ "Unknown")    
#organism <- factor(organism)

#organism <- ifelse(str_starts(age, "e") | str_starts(age, "P"), "Mouse", "Human")
#organism <- factor(organism)

#create colData with organ, age, and organism as factors
colData <- data.frame(
  organ = factor(organ),   
  age = factor(age),      
  organism = organism,           
  row.names = colnames(count_matrix_clean)
)
head(colData)
colData$age <- factor(age)
head(colData)


#check the colData
head(colData$age)

#create DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_clean, 
  colData = DataFrame(colData), 
  design = ~ organ + age
)
#apply VST normalization
dds_filtered <- dds[rowSums(counts(dds) > 0) > 0, ]
vsd <- vst(dds)

boxplot(assay(vsd)) 


pca_data <- plotPCA(vsd, intgroup = c("organ", "organism", "age"), returnData = TRUE)
#since these data are not helpful for comparing ages across organisms, I am going to
#make a new data frame for the unified ages
pca_data$age_numeric <- recode(pca_data$age,
                               "e10" = 1, "4wpc" = 1, "e93" = 1, "13" = 1, 
                               "e11" = 2, "5wpc" = 2, "e108" = 2, "14" = 2, 
                               "e12" = 3, "6wpc" = 3, "e112" = 3, "16" = 3,
                               "e13" = 4, "7wpc" = 4, "e123" = 4, "18" = 4,
                               "e14" = 5, "8wpc" = 5, "e130" = 5, "20" = 5,
                               "e15" = 6, "9wpc" = 6, "p0" = 6, "24" = 6,
                               "e16" = 7, "10wpc" = 7, "p23" = 7, "28" = 7,
                               "e17" = 8, "11wpc" = 8, "p152" = 8, "35" = 8,
                               "e18" = 9, "12wpc" = 9, "p365" = 9, "42" = 9,
                               "p0" = 10, "13wpc" = 10, "p1095" = 10, "56" = 10,
                               "p3" = 11, "16wpc" = 11, "p3285" = 11, "74" = 11,
                               "p14" = 12, "18wpc" = 12, "p5475" = 12, "104" = 12,
                               "p28" = 13, "19wpc" = 13, "p8030" = 13, "134" = 13,
                               "p63" = 14, "20wpc" = 14, "26y" = 14, "164" = 14,
                               "newborn" = 15, "p180" = 15,"194" = 15,
                               "infant" = 16,
                               "toddler" = 17, "school" = 18,
                               "youngTeenager" = 19, "youngAdult" = 20,
                               "youngMidAge" = 21, "Senior" = 22)
pca_data$age_numeric <- factor(pca_data$age_numeric)

#pca_data$age_numeric <- as.numeric(as.character(pca_data$age_numeric))
#summary(pca_data$age_numeric)
unique(pca_data$organ)
pca_data$organ_group <- recode(
  pca_data$organ, 
  "Brain" = "#1f77b4",       
  "Cerebellum" = "#1f97b4",  
  "Heart" = "#ff9f6e",      
  "Testis" = "#ff5f6e", 
  "Ovary" = "#ff5f6e",
  "Liver" = "#9ba52c",       
  "Kidney" = "#2ca02c" 
)
#pca_data$age_numeric <- as.numeric(as.character(pca_data$age_numeric))
#pca_data$age_numeric <- as.numeric(as.character(pca_data$age))
#unique(pca_data$age)
#make PCA look closer to published PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = organ_group, shape = organism, size = as.factor(age_numeric))) +
  geom_point(alpha = 0.7) + 
  labs(
    x = paste("PC1 (", round(pca_data$percentVar[1] * 100, 1), "% variance explained)", sep = ""),
    y = paste("PC2 (", round(pca_data$percentVar[2] * 100, 1), "% variance explained)", sep = ""),
    title = "Gene Expression by Organ, Organism, and Age"
  ) +
  theme_minimal() + 
  theme(
    legend.position = "right",  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12)  
  ) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +  # You can adjust this based on the number of unique organisms
  scale_size_discrete(range = c(3, 8)) +
  scale_color_identity()

