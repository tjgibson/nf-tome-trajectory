#!/usr/bin/env Rscript

# Setup ========================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(FNN)
})


options(future.globals.maxSize = 32000 * 1024^2)

# read command line arguments ==================================================
args <- commandArgs(TRUE)
time_i <- args[1]
time_j <- args[2]
embeddings_fn <- args[3]
anno_fn <- args[4]
permutation_times <-  as.numeric(args[5])
k_neigh = as.numeric(args[6])


# Read in data =================================================================
emb <- readRDS(embeddings_fn)
anno_list <- readRDS(anno_fn) |> 
  rownames_to_column("barcode") |> 
  group_by(stage) |> 
  group_split() |> 
  map(column_to_rownames, "barcode")

anno1 <- anno_list[[1]]
anno2 <- anno_list[[2]]

# Permute annotation labels for each connection ================================
res = list()

for(rep_i in 1:permutation_times){
  
  anno1$state = anno1$Anno[sample(1:nrow(anno1))]
  anno2$state = anno2$Anno[sample(1:nrow(anno2))]
  
  anno = rbind(anno1, anno2)
  if(nrow(emb) != nrow(anno)){
    stop("embeddings and annotations have different number of rows")
  }
  pd = anno[rownames(emb),]
  
  emb_sub = emb
  pd_sub = pd
  
  irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
  irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
  pd_sub1 <- pd_sub[pd_sub$day == "pre",]
  pd_sub2 <- pd_sub[pd_sub$day == "nex",]
  
  pre_state_min = min(table(as.vector(pd_sub1$state)))
  
  if (pre_state_min < k_neigh & pre_state_min >= 3){
    k_neigh = pre_state_min
    print(k_neigh)
  }
  
  if (pre_state_min < 3){
    next
  }
  
  neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
  
  tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
  for(i in 1:k_neigh){
    tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
  }
  state1 <- names(table(as.vector(pd_sub1$state)))
  state2 <- names(table(as.vector(pd_sub2$state)))
  
  tmp2 <- matrix(NA,length(state2),length(state1))
  for(i in 1:length(state2)){
    x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
    for(j in 1:length(state1)){
      tmp2[i,j] <- sum(x==state1[j])
    }
  }
  tmp2 <- tmp2/apply(tmp2,1,sum)
  tmp2 <- data.frame(tmp2)
  row.names(tmp2) = state2
  names(tmp2) = state1
  
  res[[rep_i]] = tmp2
  
}

saveRDS(res, paste0(time_i, "_", time_j, "_Knn_umap_permutation.rds"))

