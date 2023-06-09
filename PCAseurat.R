# The goal of this script is to perform PCA and clustering analysis
# on the single cell data. The tasks include:
# - Goal 1: Set the variable genes as the grin genes and cluster all cells
#    What percentage of the mixed clusters are dorsal horn neurons?
# - Goal 2: Set the variable genes to glutamate receptor genes and cluster all cells
# - Goal 3: Set the variable genes to glutamate receptor genes and cluster dorsal horn cells
# - Goal 4--: Set the variable genes as the grin genes and cluster dorsal horn neurons
# - Goal 5: Repeate Goal 1 without Grin3b or Grin2c


library(Seurat)
library(Matrix)
library(tidyverse)

# ---- Goal 1 ----
# Load single cell RNA seq data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../GScpter/Scripts/Pre_analysis_functions.R")
RDfile <- load_data("../Datasets/neuron_and_glia_2022/final_meta_dataset.rds")

# Remove the counts to start a fresh with a new Seurat object
normed_data <- GetAssayData(RDfile, assay = "raw", slot = "data")
meta_data <- RDfile@meta.data
rm(RDfile)
gc()


# Set the variable features
varfeatures <- c("Grin1", "Grin2a", "Grin2b", "Grin2c", "Grin2d", "Grin3a", "Grin3b")
# Now run through the clustering pipeline
# Skip the QC steps because they typically have already been done
# check if the counts have decimals and see if it needs normalizing
# Just take the data one and use it from the normalization step
# and also what normalization are we doing and just start with scaling.
(neurnastro <- which(meta_data$final_coarse_types %in% c("Neurons", "Astrocytes")))
length(neurnastro)
grinclu <- CreateSeuratObject(counts = normed_data[, neurnastro],
                              meta.data = meta_data[neurnastro, ])
(all.genes <- rownames(grinclu))
grinclu <- ScaleData(grinclu, all.genes)
grinclu <- FindVariableFeatures(grinclu)
# Which of the Grin genes is a variable feature of the astrocytes and neurons
varfeatures[which(varfeatures %in% grinclu@assays$RNA@var.features)]
# Only Grin2a.
# What about for all ionotropic receptors?
var2features <- c("Grin1", "Grin2a", "Grin2b", "Grin2c", "Grin2d", "Grin3a", "Grin3b",
                 "Grik1", "Grik2", "Grik3", "Grik4", "Grik5",
                 "Gria1", "Gria2", "Gria3", "Gria4")
(vars <- var2features[which(var2features %in% grinclu@assays$RNA@var.features)])
# Only 5 of the ionotropic glutamate receptors are considered variable features

grinclu <- RunPCA(grinclu)
grinclu <- RunPCA(grinclu, features = varfeatures)
grinclu <- RunPCA(grinclu, features = var2features)
grinclu <- RunPCA(grinclu, features = vars)
# I get a bunch of warnings as I thought
DimPlot(grinclu, reduction = "pca", group.by = "final_coarse_types",
        pt.size = 1)
DimPlot(grinclu, reduction = "pca", group.by = "final_cluster_assignment")
ElbowPlot(grinclu)

# Clustering
grinclu <- FindNeighbors(grinclu, dims = 1:6)
grinclu <- FindClusters(grinclu, resolution = 2)
head(Idents(grinclu), 5)

grinclu <- RunUMAP(grinclu, dims = 1:6)
DimPlot(grinclu, reduction = "umap")
DimPlot(grinclu, reduction = "umap", group.by = "final_coarse_types")
DimPlot(grinclu, reduction = "umap", group.by = "final_cluster_assignment")
DotPlot(grinclu, features = varfeatures)
DotPlot(grinclu, features = varfeatures, group.by = "final_coarse_types")
DotPlot(grinclu, features = varfeatures, group.by = "final_cluster_assignment")
FeaturePlot(grinclu, features = varfeatures)

saveRDS(grinclu, "goal1.rds")
end <- Sys.time()
timing <- c(timing, goal1 = end - start)
goal1 <- readRDS("goal1.rds")
(pca <- DimPlot(goal1, reduction = "pca"))
ggsave("goal1pca.png", plot = pca, height = 1200, device = "png",
       width = 2100, units = "px", dpi = 300, type = "cairo")
(umap <- DimPlot(goal1, reduction = "umap"))
ggsave("goal1umap.png", plot = umap, height = 1200, device = "png",
       width = 2100, units = "px", dpi = 300, type = "cairo")
(dot <- DotPlot(goal1, features = varfeatures))
ggsave("goal1dotplot.png", plot = dot, height = 1200, device = "png",
       width = 2100, units = "px", dpi = 300, type = "cairo")

# ---- Goal 2 ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../GScpter/Scripts/Pre_analysis_functions.R")
RDfile <- load_data("../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")

# Remove the counts to start a fresh with a new Seurat object
start <- Sys.time()
normed_data <- GetAssayData(RDfile, assay = "raw", slot = "data")
rm(RDfile)
gc()

# Set the variable features
varfeatures <- c("Grin1", "Grin2a", "Grin2b", "Grin2c", "Grin2d", "Grin3a", "Grin3b",
                 "Grik1", "Grik2", "Grik3", "Grik4", "Grik5",
                 "Gria1", "Gria2", "Gria3", "Gria4")

# Now run through the clustering pipeline
# Skip the QC steps because they typically have already been done
# check if the counts have decimals and see if it needs normalizing
# Just take the data one and use it from the normalization step and
# and also what normalization are we doing and just start with scaling.
grinclu <- CreateSeuratObject(counts = normed_data)
rm(normed_data)
gc()
all.genes <- rownames(grinclu)
grinclu <- ScaleData(grinclu, all.genes)
grinclu <- RunPCA(grinclu, features = varfeatures)
# I get a bunch of warnings as I thought
DimPlot(grinclu, reduction = "pca")
ElbowPlot(grinclu)

# Clustering
grinclu <- FindNeighbors(grinclu, dims = 1:15)
grinclu <- FindClusters(grinclu, resolution = 0.5)
head(Idents(grinclu), 15)

grinclu <- RunUMAP(grinclu, dims = 1:6)
DimPlot(grinclu, reduction = "umap")
saveRDS(grinclu, "goal2.rds")
end <- Sys.time()
timing <- c(timing, goal2 = end - start)
goal2 <- readRDS("goal2.rds")
(dot <- DotPlot(goal2, features = varfeatures))


# ---- Goal 3 ----

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../GScpter/Scripts/Pre_analysis_functions.R")
RDfile <- load_data("../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")

# Get the dorsal horn neurons
start <- Sys.time()
all_clusters <- c("Excit-1", "Excit-2", "Excit-3", "Excit-4", "Excit-5",
                  "Excit-6", "Excit-8", "Excit-9","Excit-10","Excit-12",
                  "Excit-13","Excit-14","Excit-15","Excit-16","Excit-18",
                  "Excit-19","Excit-20","Inhib-1","Inhib-2","Inhib-3",
                  "Inhib-4","Inhib-5","Inhib-6","Inhib-7","Inhib-8","Inhib-9",
                  "Inhib-10","Inhib-11","Inhib-12","Inhib-13","Inhib-14",
                  "Excit-21","Excit-22","Excit-23","Excit-24","Excit-25",
                  "Excit-26","Excit-27","Excit-29","Excit-30","Inhib-15",
                  "Inhib-16","Inhib-18","Inhib-19","Inhib-20","Inhib-21",
                  "Excit-31","Excit-32","Excit-34","Excit-35","Excit-36")
CellIndicies <- which(meta_ident[["final_clusters"]] %in% unique(all_clusters))

# Remove the counts to start a fresh with a new Seurat object
normed_data <- GetAssayData(RDfile, assay = "raw", slot = "data")[, CellIndicies]
rm(RDfile)
gc()

# Set the variable features
varfeatures <- c("Grin1", "Grin2a", "Grin2b", "Grin2c", "Grin2d", "Grin3a", "Grin3b",
                 "Grik1", "Grik2", "Grik3", "Grik4", "Grik5",
                 "Gria1", "Gria2", "Gria3", "Gria4")

# Now run through the clustering pipeline
# Skip the QC steps because they typically have already been done
# check if the counts have decimals and see if it needs normalizing
# Just take the data one and use it from the normalization step and
# and also what normalization are we doing and just start with scaling.
grinclu <- CreateSeuratObject(counts = normed_data)
rm(normed_data)
gc()
all.genes <- rownames(grinclu)
grinclu <- ScaleData(grinclu, all.genes)
grinclu <- RunPCA(grinclu, features = varfeatures)
# I get a bunch of warnings as I thought
DimPlot(grinclu, reduction = "pca")
ElbowPlot(grinclu)

# Clustering
grinclu <- FindNeighbors(grinclu, dims = 1:6)
grinclu <- FindClusters(grinclu, resolution = 0.5)
head(Idents(grinclu), 5)

grinclu <- RunUMAP(grinclu, dims = 1:6)
DimPlot(grinclu, reduction = "umap")
saveRDS(grinclu, "goal3.rds")
end <- Sys.time()
timing <- c(timing, goal3 = end - start)



# ---- Goal 4 ----

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../GScpter/Scripts/Pre_analysis_functions.R")
RDfile <- load_data("../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")

# Get the dorsal horn neurons

all_clusters <- c("Excit-1", "Excit-2", "Excit-3", "Excit-4", "Excit-5",
                  "Excit-6", "Excit-8", "Excit-9","Excit-10","Excit-12",
                  "Excit-13","Excit-14","Excit-15","Excit-16","Excit-18",
                  "Excit-19","Excit-20","Inhib-1","Inhib-2","Inhib-3",
                  "Inhib-4","Inhib-5","Inhib-6","Inhib-7","Inhib-8","Inhib-9",
                  "Inhib-10","Inhib-11","Inhib-12","Inhib-13","Inhib-14",
                  "Excit-21","Excit-22","Excit-23","Excit-24","Excit-25",
                  "Excit-26","Excit-27","Excit-29","Excit-30","Inhib-15",
                  "Inhib-16","Inhib-18","Inhib-19","Inhib-20","Inhib-21",
                  "Excit-31","Excit-32","Excit-34","Excit-35","Excit-36")
CellIndicies <- which(meta_ident[["final_clusters"]] %in% unique(all_clusters))

# Remove the counts to start a fresh with a new Seurat object
normed_data <- GetAssayData(RDfile, assay = "raw", slot = "data")[, CellIndicies]
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
rm(normed_data)
gc()
all.genes <- rownames(grinclu)
grinclu <- ScaleData(grinclu, all.genes)
grinclu <- RunPCA(grinclu, features = varfeatures)
# I get a bunch of warnings as I thought
DimPlot(grinclu, reduction = "pca")
ElbowPlot(grinclu)

# Clustering
grinclu <- FindNeighbors(grinclu, dims = 1:6)
grinclu <- FindClusters(grinclu, resolution = 0.5)
head(Idents(grinclu), 5)

grinclu <- RunUMAP(grinclu, dims = 1:6)
DimPlot(grinclu, reduction = "umap")
saveRDS(grinclu, "goal4.rds")
end <- Sys.time()
timing <- c(timing, goal4 = end - start)



# ---- Goal 5a ----
# cluster all cells with out grin 2c and grin3b
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../GScpter/Scripts/Pre_analysis_functions.R")
RDfile <- load_data("../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")

# Remove the counts to start a fresh with a new Seurat object
normed_data <- GetAssayData(RDfile, assay = "raw", slot = "data")
rm(RDfile)
gc()

# Set the variable features
varfeatures <- c("Grin1", "Grin2a", "Grin2b", "Grin2d", "Grin3a")

# Now run through the clustering pipeline
# Skip the QC steps because they typically have already been done
# check if the counts have decimals and see if it needs normalizing
# Just take the data one and use it from the normalization step and
# and also what normalization are we doing and just start with scaling.
grinclu <- CreateSeuratObject(counts = normed_data)
rm(normed_data)
gc()
all.genes <- rownames(grinclu)
grinclu <- ScaleData(grinclu, all.genes)
grinclu <- RunPCA(grinclu, features = varfeatures)
# I get a bunch of warnings as I thought
DimPlot(grinclu, reduction = "pca")
ElbowPlot(grinclu)

# Clustering
grinclu <- FindNeighbors(grinclu, dims = 1:6)
grinclu <- FindClusters(grinclu, resolution = 0.5)
head(Idents(grinclu), 5)

grinclu <- RunUMAP(grinclu, dims = 1:6)
DimPlot(grinclu, reduction = "umap")
saveRDS(grinclu, "goal5a.rds")



# ---- Goal 5b ----
# remove the genes with least variance
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../GScpter/Scripts/Pre_analysis_functions.R")
RDfile <- load_data("../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")

# Get the dorsal horn neurons
all_clusters <- c("Excit-1", "Excit-2", "Excit-3", "Excit-4", "Excit-5",
                  "Excit-6", "Excit-8", "Excit-9","Excit-10","Excit-12",
                  "Excit-13","Excit-14","Excit-15","Excit-16","Excit-18",
                  "Excit-19","Excit-20","Inhib-1","Inhib-2","Inhib-3",
                  "Inhib-4","Inhib-5","Inhib-6","Inhib-7","Inhib-8","Inhib-9",
                  "Inhib-10","Inhib-11","Inhib-12","Inhib-13","Inhib-14",
                  "Excit-21","Excit-22","Excit-23","Excit-24","Excit-25",
                  "Excit-26","Excit-27","Excit-29","Excit-30","Inhib-15",
                  "Inhib-16","Inhib-18","Inhib-19","Inhib-20","Inhib-21",
                  "Excit-31","Excit-32","Excit-34","Excit-35","Excit-36")
CellIndicies <- which(meta_ident[["final_clusters"]] %in% unique(all_clusters))

# Remove the counts to start a fresh with a new Seurat object
normed_data <- GetAssayData(RDfile, assay = "raw", slot = "data")[, CellIndicies]
rm(RDfile)
gc()

# Set the variable features
varfeatures <- c("Grin1", "Grin2a", "Grin2b", "Grin2d", "Grin3a", "Grin3b",
                 "Grik1", "Grik2", "Grik3", "Grik4", "Grik5",
                 "Gria1", "Gria2", "Gria3", "Gria4")

# Now run through the clustering pipeline
# Skip the QC steps because they typically have already been done
# check if the counts have decimals and see if it needs normalizing
# Just take the data one and use it from the normalization step and
# and also what normalization are we doing and just start with scaling.
grinclu <- CreateSeuratObject(counts = normed_data)
rm(normed_data)
gc()
all.genes <- rownames(grinclu)
grinclu <- ScaleData(grinclu, all.genes)
grinclu <- RunPCA(grinclu, features = varfeatures)
# I get a bunch of warnings as I thought
DimPlot(grinclu, reduction = "pca")
ElbowPlot(grinclu)

# Clustering
grinclu <- FindNeighbors(grinclu, dims = 1:6)
grinclu <- FindClusters(grinclu, resolution = 0.5)
head(Idents(grinclu), 5)

grinclu <- RunUMAP(grinclu, dims = 1:6)
DimPlot(grinclu, reduction = "umap")
saveRDS(grinclu, "goal5b.rds")





# ---- DE and trends ----
# Differential gene expression and trends
# Make another script for this
