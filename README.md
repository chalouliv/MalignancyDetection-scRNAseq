# MalignancyDetection-scRNAseq
COMP402: Year Long Honors Biology and Computer Science Project. 

## Testing a Novel tool for Malignancy Determination in Single-Cell RNA-Seq in Pediatric Brain Cancer
Author: Charlotte Livingston 
Principal Investigator: Claudia Kleinman 
Date: April 27th, 2025

### Motivation
Tumors are complex where malignant and non-malignant cells interact within the tumor microenvironment (TME). These interactions influence key processes such as angiogenesis, proliferation, metastasis, and treatment resistance. Single-cell RNA sequencing (scRNA-seq) has transformed our understanding of the TME by revealing cellular heterogeneity at single-cell resolution, enabling improved tumor classification and progression prediction.

A major challenge in analyzing scRNA-seq tumor data is distinguishing malignant from non-malignant cells based on gene expression. This task is complicated by low sequencing coverage and transcriptional similarity between cell types. This is particularly true in pediatric brain tumors, where malignant cells often resemble their non-malignant counterparts.

Copy number variation (CNV) inference from gene expression is the most common method for identifying malignancy, but it depends on accurate references and manual interpretation, leading to potential bias and reduced sensitivity for rare cell populations.

Consensus Reference-based Automated Labeling (CoRAL), developed by the Kleinman lab, integrates over 20 machine learning and statistical annotation methods to produce consensus-based cell-type labels. My COMP402 project investigates whether CoRAL can also accurately distinguish malignant from normal cells in scRNA-seq data, providing a more robust and automated approach to malignancy detection.

### Methods
CoRAL's effectiveness was evaluated in human pediatric tumors aswell as mice allografts.
