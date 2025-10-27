# ------------------------------------------------------------
# Script Name: combine_human.R
# Description: Combining seurat objects to make human training set 
# Author: Charlotte Livingston
# Email: charlotte.livingston@mail.mcgill.ca
# Date: 2024-12-01
# Project: COMP402 Honors Project 
# ------------------------------------------------------------

setwd("/path/")
source("functions.R")

## Load objects 

# Updating and loading etmr1 object
ETMR1 <- loadRDa("/path/") 
ETMR1@reductions$pca@global <- TRUE
ETMR1@reductions$tsne@global <- TRUE
ETMR1 <- UpdateSeuratObject(ETMR1)

# Loading other tumors
BT2016062_PFA <- loadRDa("/path/") 
CNS_NB_GG <- loadRDa("/path/") 
HSJ_152sc <- loadRDa("/path/") 

BT2016062_PFA <- UpdateSeuratObject(BT2016062_PFA)
CNS_NB_GG <- UpdateSeuratObject(CNS_NB_GG)
HSJ_152sc <- UpdateSeuratObject(HSJ_152sc)

# Renaming metadata columns
BT2016062_PFA$cell_type <- BT2016062_PFA$old_column
HSJ_152sc$cell_type <- HSJ_152sc$old_column

# Loading healthy samples
kim <- loadRDa("/path/") 
lin <- loadRDa("/path/")
vel <- loadRDa("/path/") 

# Combine all objects 
set_1 <- merge(BT2016062_PFA, y = c(ETMR1, CNS_NB_GG, HSJ_152sc, kim, lin, vel), add.cell.ids = c("A", "B", "C", "D", "E", "F", "G"))  %>% add_malignancy() %>% add_broad_celltype() 

# Add metadata
set_1 <- add_normal_celltypes(set_1) 
set_1 <- add_normal_celltypes(set_1) 

# Save
save(set_1, file = "/path/") 

# Removing 'Uncertain' Cells for clean benchmarking 
set_1_wo <- set_1[,(set_1$normal_celltype %in% c('Malignant', 'Normal.Glia', 'Normal.Immune', 'Normal.Neurons', 'Normal.Vascular'))]
save(set_1_wo, file = "/path/") 

