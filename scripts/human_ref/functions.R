## function for loading .rda ##
loadRDa <- function(fileName){
  # Loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## function for loading sample list ##
get_samples <- function(samples) {
  # Use lapply to retrieve each sample by name and return as a list
  return(lapply(samples, get))
}

# ------------------------------------
### FUNCTIONS TO HARMONIZE NAMING ###
# ------------------------------------

# 1 of 2 functions for harmonizing celltype label
board_celltype_labeling <- function(celltype){
  # takes an old celltype and returns the harmonized label
  if(celltype %in% c("Vascular & other")) {
    return("Vascular")
  }
  else if (celltype %in% c("Neurons")) {
    return("Neurons")
  }
  else if (celltype %in% c("Astrocytes", "Ependymal", "Oligodendrocytes", "OPC", "Proliferating OPC")) {
    return("Glia")
  }
  else if (celltype %in% c("Immune")) {
    return("Immune")
  }
  else{
    return("Uncertain")
  }
}

## 2 of 2 functions for harmonizing celltype label per cell across objects ##
add_broad_celltype <- function(sample){
  # Creates new meta_data column "broad_celltype" with harmonized celltype labels
  sample@meta.data$broad_celltype <- mapply(function(class, consensus) {
    if (!is.na(class)){
      class
    }
    else if(!is.na(consensus)){
      board_celltype_labeling(consensus)
    }
    else {
      "Uncertain"
    }
  }, sample$class, sample$cell_type)
  return(sample)
}

# Function to add combined Malignancy + Celltype labels to only normal cells 
add_normal_celltypes <- function(sample){
  # Creates new meta_data column "normal_celltype" with harmonized 'Normal.celltype' label 
  # for normal cells and 'Malignant' for malignant cells
  sample@meta.data$normal_celltype <- mapply(function(malignancy, broad_celltype) {
    if (malignancy == "Normal") { 
      paste0("Normal.", broad_celltype)
    }
    else{
      malignancy
    }
  }, sample$Malignancy, sample$broad_celltype)
  return(sample)
}

# Function to add combined Malignancy + Celltype labels to only normal cells  
add_all_celltypes <- function(sample){
  # Creates new meta_data column "malignant_celltype" with harmonized label for all cells
  sample@meta.data$malignancy_celltype <- mapply(function(malignancy, celltype) {
    paste0(malignancy, ".", celltype)
  }, sample$malignancy, sample$broad_celltype)
  return(sample)
}

# Function to evaluate prediction correctness per cell
correctness <- function(sample, true_col) {
  # Creates new meta_data column "orrect_Predicition" with...
  # "Agree" if the prediction exactly matches the ground truth 
  # or "Disagree" elsewise
  sample@meta.data$Correct_Predicition <- mapply(function(true, pred) {
    if (true == pred) {
      return("Agree")
    }
    else {
      return("Disagree")
    }
  }, sample@meta.data[[truth_col]], sample$CAWPE_Prediction)
  return(sample)
}

# ------------------------------------
### FUNCTIONS FOR PLOTTING ###
# ------------------------------------

# Plot confusion matrix as a heatmap 
plot_cm = function(cm_table){
  col_fun = circlize::colorRamp2(c(range(cm_table)[1], 
                                   range(cm_table)[2]/2, 
                                   range(cm_table)[2]), 
                                 c("#5C80BC", "#F2EFC7", "#FF595E")) 
  
  h = Heatmap(cm_table,
              name = 'Counts',
              col = col_fun,
              width = ncol(cm_table)*unit(6, "mm"),
              height = nrow(cm_table)*unit(6, "mm"),
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 15), 
              column_title = 'True Class', 
              row_title = 'Predicted Class')
  
  return(h)
}

custom_v <- function(sample, label, feature) {
  
  # Create violin
  v <- VlnPlot(sample, features = feature, group.by = label, pt.size = 0) +
    theme(aspect.ratio = 1, legend.position = "none", axis.text.x = element_text(size = 8)) +
    labs(x = "T: Ground Truth Label \n P: Predicted Label", y = "Expression Level")
  
  # to add cell counts
  data <- v$data
  colnames(data) <- c("exp", "ident")
  counts <- as.data.frame(table(data$ident))
  colnames(counts) <- c("ident", "Freq")
  
  # Get height per group
  y_max <- aggregate(data$exp, by = list(data$ident), FUN = max)
  colnames(y_max) <- c("ident", "y_max")
  
  # Merge to prepare for labeling
  label_data <- merge(counts, y_max, by = "ident")
  
  # Create a geom_text layer for cell counts
  c <- geom_label(data = label_data,
                  aes(x = ident, y = y_max * 1.1, label = Freq),
                  fill = "white",     # White background
                  color = "black",    # Black text color
                  size = 4,
                  label.size = 0)
  
  #extend y-axis to include space for the counts label 
  max <- max(y_max$y_max)
  
  return(v+c+ylim(0, max * 1.2))
}