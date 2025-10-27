# function for loading .rda ----------------
loadRDa <- function(fileName){
  # Loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

get_samples <- function(samples) {
  # Use lapply to retrieve each sample by name and return as a list
  return(lapply(samples, get))
}

# ------------------------------------
### FUNCTIONS TO HARMONIZE NAMING ###
# ------------------------------------

# Function to add malignancy labeling 
add_malignancy <- function(sample, malignant, normal, uncertain){
  sample@meta.data$Malignancy <- mapply(function(infer_call, gfp, doublet) {
    # Define new column (malignancy) based on inferCNV and GFP
    if (doublet == "doublet") {
      "Doublet"
    }
    else if (is.element(infer_call, normal)) { # Normal inferCNVs
      if(gfp == "negative") {
        "Normal"}
      else {
        "Uncertain"}     
    }
    else if(is.element(infer_call, malignant))  { #malignant inferCNV
      if(gfp == "positive") {
        "Malignant"}
      else {
        "Uncertain"}
    }
    else {
      "Uncertain"}
  }, sample$inferCNV, sample$GFP_status, sample$doublet_call)
  return(sample)
}

# change GFP postive v. negative to follow less strict requirements
fix_gfp <- function(sample){
  gfp_expression <- Seurat::FetchData(sample, vars = "GFP")
  gfp_expression_class <-ifelse(gfp_expression > 0.48, "positive", "negative")
  sample@meta.data[["GFP_status"]] <- gfp_expression_class
  return(sample)
}

# Function to add celltype labeling  to normal cells 
add_normal_celltypes <- function(sample){
  sample@meta.data$normal_celltype <- mapply(function(malignancy, celltype) {
    # Define new column (normal_celltype) based on inferCNV and GFP
    if (malignancy == "normal") {
      paste0("Normal.",celltype)
    }
    else{
      "Malignant"
    }
  }, sample$malignant, sample$broad_celltype)
  return(sample)
}

# Function to add celltype labeling  to all cells 
add_all_celltypes <- function(sample){
  sample@meta.data$malignancy_celltype <- mapply(function(malignancy, celltype) {
    # Define new column (malignancy_celltype) based on inferCNV and GFP
    paste0(malignancy, ".", celltype)
  }, sample$Malignancy, sample$broad_celltype)
  return(sample)
}

# Function to add broad celltype labeling
broad_celltype_labeling <- function(celltype) {
  if (celltype %in% c("Vascular & other", "Vascular", "Pericytes", "Endothelial")) {
    return("Vascular_other")
  }
  else if (celltype %in% c("Neurons", "Neuronal_progenitors", "CGE Inhibitory", "Excitatory L2-L4", "Excitatory L5-L6", "MGE Inhibitory", "Other Excitatory", "RGC", "Other Inhibitory")) {
    return("Neurons")
  }
  #else if (celltype %in% c("Astrocytes", "Ependymal", "Oligodendrocytes", "OPC", "Proliferating OPC", "Choroid_plexus", "Glial_progenitors")) {
  #  return("Glia")
  #}
  else if (celltype %in% c("Immune", "T_Cells", "B_Cells", "Macrophages", "DC", "Monocytes", "NK_Cells", "Mast_Cells", "Microglia", "Meninges")) {
    return("Immune")
  }
  else {
    return(celltype)
  }
}

# Function to add broad celltype labels to the sample
add_broad_celltype <- function(sample) {
  sample@meta.data$broad_celltype <- mapply(function(celltype) {
    if (!is.na(celltype)) {
      broad_celltype_labeling(celltype)
    }
    else {
      return("Uncertain")
    }
  }, sample$celltype)
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


# custom violin. Takes as input seurat object, label = name of meta data column name to group.by, and feature
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