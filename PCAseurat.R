# The goal of this script is to perform PCA and clustering analysis
# on the single cell data. The tasks include:
# - Goal 1: Set the variable genes as the grin genes and cluster all cells
# - Goal 2: Set the variable genes as the grin genes and cluster dorsal horn cells

library(Seurat)
library(tidyverse)

# Goal 1
# Load single cell RNA seq data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../GScpter/Scripts/Pre_analysis_functions.R")
RDfile <- load_data("../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")

# Remove the counts to start a fresh with a new Seurat object
normed_data <- GetAssayData(RDfile, assay = "RNA", slot = "data")
rm(RDfile)
gc()

# Set the variable features
varfeatures <- c("Grin1", "Grin2a", "Grin2b", "Grin2c", "Grin2d", "Grin3a", "Grin3b")

# Now run through the clustering pipeline
# Skip the QC steps because they typically have already been done
# check if the counts have decimals and see if it needs normalizing
# Just take the data one and use it from the normalization step and
# and also what normalization are we doing and just start with scaling.
grinclu <- CreateSeuratObject(counts = normed_data)
all.genes <- rownames(grinclu)
grinclu <- ScaleData(grinclu, all.genes)
grinclu <- RunPCA(grinclu, features = varfeatures)