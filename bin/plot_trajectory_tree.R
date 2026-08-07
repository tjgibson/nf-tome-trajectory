#!/usr/bin/env Rscript

# Setup ========================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(tidygraph)
  library(ggraph)
})

# read command line arguments ==================================================
args <- commandArgs(TRUE)
tp_list <- unlist(strsplit(args[1], split = ","))
edge_prob_thresh <- as.numeric(args[2])
edge_all_fn <- args[3]

# read in edge list ============================================================
edge_all <- readRDS(edge_all_fn)

# convert to format for ggraph =================================================
edge_data <-  edge_all |> 
  dplyr::filter(prob > edge_prob_thresh) |> 
  dplyr::select(1:3) |> 
  dplyr::rename(from = pre, to = nex, weight = prob) |> 
  dplyr::select(from, to, weight) 

node_data <- tibble(
  node_name = unique(c(edge_data$from, edge_data$to)) 
) |> 
  tidyr::separate_wider_delim(node_name, delim = ":", names = c("stage", "cell_type"), cols_remove = FALSE)

p_data <- tbl_graph(
  edges = edge_data, 
  nodes = node_data,
  node_key = "node_name"
) 

# plot tree ====================================================================
n_tp <- length(tp_list)
n_celltype <- length(unique(node_data$cell_type))

pdf("tree.pdf", useDingbats = FALSE, width = n_celltype * 1.5, height = n_tp * 1.5)

ggraph(p_data, "tree") +
  geom_edge_diagonal(aes(alpha = weight)) +
  geom_node_point(aes(color = cell_type)) +
  geom_node_text(aes(label = node_name))

dev.off()
