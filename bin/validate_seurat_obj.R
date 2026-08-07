#!/usr/bin/env Rscript

# Setup ========================================================================
suppressPackageStartupMessages({
  library(Seurat)
})

# read command line arguments ==================================================
args <- commandArgs(TRUE)
seurat_object_fn <- args[1]
timepoint_metadata_key <- args[2]
tp_list <- unlist(strsplit(args[3], split = ","))

# Read in data =================================================================
seurat_obj <- readRDS(seurat_object_fn)
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj <- UpdateSeuratObject(seurat_obj)

# perform input validation =====================================================
if (!(timepoint_metadata_key %in% colnames(seurat_obj@meta.data))) {
  stop(paste("timepoint_metadata_key", timepoint_metadata_key, "is not a column in the metadata for input seurat object"))
}

if (any(!(tp_list %in% seurat_obj@meta.data[[timepoint_metadata_key]]))) {
  stop("not all timepoint_names are present in the seurat object")
}

if (!("RNA" %in% Assays(seurat_obj))) {
  stop("Seurat object is missing RNA assay")
}
  
  
  