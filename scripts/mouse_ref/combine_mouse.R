# ------------------------------------------------------------
# Script Name: combine_mouse.R
# Description: Combining seurat objects to make human training set 
# Author: Charlotte Livingston
# Email: charlotte.livingston@mail.mcgill.ca
# Date: 2024-12-01
# Project: COMP402 Honors Project 
# ------------------------------------------------------------

setwd("/path/")
source("functions.R")

## Load objects 
# Tumor objects
samples <- c("7944_GPAC", "8008_KPP", "AN14895_KPF", "25097_KP", "19074_KPPMPIK", "25091_KNF", "7961_KF", "5-1g_GPAC", "7992_EZHIP_p53", "16-1A_H31KP", "8007_KPP", "13-1a_EZHIP_p53", "3-3g_H31KP")
s_path <-  "/path/"
for (s in samples){
  temp <- loadRDa(paste0(s_path, s, ".Rda"))
  temp <- UpdateSeuratObject(temp)
  temp <- add_broad_celltype(temp) # Harmonize celltype labelling
  temp <- add_malignancy(temp) # Harmonize malignancy labelling
  assign(s, temp, envir = .GlobalEnv)
}

# Adding suplementary normal cells
normal_cells <- loadRDa("/path/")

# Randomly sample cells from the normal object
random_cells <- sample(Cells(normal_cells), 10000)
normal_cells <- subset(normal_cells, cells = random_cells)

# Compare features to ensure comparable quality 
Plot1 <- Seurat::VlnPlot(set_2, c('nCount_RNA', 'nFeature_RNA'), ncol = 2, pt.size = -1, group.by = 'orig.ident')
ggplot2::ggsave("/path/", plot = Plot2, width = 6, height = 4, dpi = 300)

Plot2 <- Seurat::VlnPlot(normal_cells, c('nCount_RNA', 'nFeature_RNA'), ncol = 2, pt.size = -1, group.by = 'orig.ident')
ggplot2::ggsave("/path/", plot = Plot2, width = 6, height = 4, dpi = 300)

# Add meta_data to normal cells
normal_cells@meta.data$Malignancy <- "Normal"
normal_cells <- add_celltype_normal(normal_cells)

# Combine normal and tumor samples
set_2 <- normal_cells
set_2 <- merge(x= set_2, y = get_samples(samples), add.cell.ids = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"))

# REMOVE GFP from expression matrix (dirty)
set_2[["RNA"]]@counts["GFP", ] <- 0
set_2[["RNA"]]@data["GFP", ] <- 0

# Add the metadata
set_2 <- add_all_celltypes(set_2)
set_2 <- add_normal_celltypes(set_2)
list <- unique(set_2$normal_celltype)
set_2 <- set_2[,(set_2$normal_celltype %in% list)]

save(set_2, file = "/path/") 
