#!/usr/bin/env Rscript

# Setup ========================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})


options(future.globals.maxSize = 32000 * 1024^2)

# read command line arguments ==================================================
args <- commandArgs(TRUE)
seurat_object_fn <- args[1]
timepoint_metadata_key <- args[2]
cell_type_key <- args[3]
time_1 <- args[4]
time_2 <- args[5]
lineage_knn <- args[6]
k_neighbors <- as.numeric(args[7])

source(lineage_knn)
# Read in data =================================================================
message("reading seurat object")
seurat_obj <- readRDS(seurat_object_fn)
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj <- UpdateSeuratObject(seurat_obj)

# subset to timepoints of interest =============================================
tp_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[timepoint_metadata_key]] %in% c(time_1, time_2),])
seurat_tp_subset <- subset(seurat_obj, cells = tp_cells)
rm(seurat_obj)

seurat_tp_subset[["RNA"]] <- split(seurat_tp_subset[["RNA"]], f = seurat_tp_subset@meta.data[[timepoint_metadata_key]])

# clustering and integration of selected timepoints ============================
message("performing clustering and integration")
seurat_tp_subset <- NormalizeData(seurat_tp_subset, verbose = FALSE)
seurat_tp_subset <- FindVariableFeatures(seurat_tp_subset, verbose = FALSE)
seurat_tp_subset <- ScaleData(seurat_tp_subset, verbose = FALSE)
seurat_tp_subset <- RunPCA(seurat_tp_subset, verbose = FALSE)

seurat_tp_subset <- IntegrateLayers(
  object = seurat_tp_subset, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  verbose = FALSE
)

seurat_tp_subset <- ScaleData(seurat_tp_subset, verbose = FALSE)
# seurat_tp_subset <- FindNeighbors(seurat_tp_subset, reduction = "integrated.rpca", dims = 1:cluster_npcs)
# seurat_tp_subset <- FindClusters(seurat_tp_subset, resolution = 1, cluster.name = "rpca_clusters")

seurat_tp_subset <- RunUMAP(
  seurat_tp_subset, 
  reduction = "integrated.rpca", 
  dims = 1:30, 
  n.components = 3, 
  min.dist = 0.75,
  reduction.name = "umap.rpca"
  )

# write seurat object for pairwise integration =================================
message("saving results pairwise integration")
saveRDS(seurat_tp_subset, paste0(time_1, "_", time_2, ".rds"))

# write embeddings to file =====================================================
emb <-  data.frame(Embeddings(object = seurat_tp_subset, reduction = "umap.rpca"))
saveRDS(emb, file=paste0(time_1, "_", time_2, "_umap3.rds"))

# Run Knn to find ancestor =====================================================
message("preparing annotations for createLineage_Knn")
library(FNN)
# prep data frame with stage and cell type
anno <- seurat_tp_subset@meta.data |> 
  as.data.frame() |> 
  mutate(day = case_when(
    !!sym(timepoint_metadata_key) == time_1 ~ "pre",
    !!sym(timepoint_metadata_key) == time_2 ~ "nex"
  )) |> 
  rename(stage = !!sym(timepoint_metadata_key)) |> 
  rename(Anno = !!sym(cell_type_key)) |> 
  select(day, Anno, stage)

# save to file
saveRDS(anno, paste0(time_1, "_", time_2, "_anno.rds"))

# run knn
message("running createLineageKnn")
res <-  createLineage_Knn(emb, anno,  k_neigh = k_neighbors)

# save to file 
message("saving results to file")
saveRDS(res, paste0(time_1, "_", time_2, "_Knn_umap.rds"))
