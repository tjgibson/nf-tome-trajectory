#!/usr/bin/env Rscript

# Setup ========================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(reshape2)
})


options(future.globals.maxSize = 32000 * 1024^2)

# read command line arguments ==================================================
args <- commandArgs(TRUE)
tp_list <- unlist(strsplit(args[1], split = ","))
replication_times <-  as.numeric(args[2])
edge_prob_thresh <- as.numeric(args[3])


# extract data for all edges ===================================================
res_median_umap = list()
for(time_i in 1:(length(tp_list)-1)){
  message(paste0(tp_list[time_i], ":", tp_list[time_i+1]))
  dat = readRDS(paste0(tp_list[time_i],"_",tp_list[time_i+1],"_Knn_umap.rds"))
  state_1 <- paste0(tp_list[time_i+1],":",row.names(dat[[1]]))
  state_2 = paste0(tp_list[time_i],":",colnames(dat[[1]])) 

  tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))
  for(i in 1:nrow(dat[[1]])){
    for(j in 1:ncol(dat[[1]])){
      xx = NULL
      for(k in 1:replication_times){
        xx = c(xx, dat[[k]][i,j])
      }
      tmp_1[i,j] = median(xx[!is.na(xx)])
    }
  }
  tmp_1 = data.frame(tmp_1)
  row.names(tmp_1) = state_1
  colnames(tmp_1) = state_2
  res_median_umap[[time_i]] = tmp_1
}

dat = NULL
for(i in 1:length(res_median_umap)){
  print(tp_list[i])
  dat = rbind(dat, melt(as.matrix(res_median_umap[[i]])))
}

dat = data.frame(dat)
names(dat) = c("nex", "pre", "prob")

dat$pre_time = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][1]))

# add root node
t1_nodes <- dat |> 
  filter(pre_time == tp_list[1]) |> 
  pull(pre) |> 
  as.character() |> 
  unique() 

root_edges <- tibble(
  nex = t1_nodes,
  pre = "root:root",
  prob = 1,
  pre_time = tp_list[1]
  
)

dat <- bind_rows(root_edges, dat)



dat$pre_cell = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat$nex_time = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][1]))
dat$nex_cell = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][2]))



saveRDS(dat, paste0("edge_all.rds"))


print(paste0("how many edges: ", nrow(dat)))
print(paste0("how many edges (> 0): ", nrow(dat[dat$prob>0,])))
print(paste0("how many edges (> 0.2): ", nrow(dat[dat$prob>=0.2,])))
print(paste0("how many edges (> 0.7): ", nrow(dat[dat$prob>=0.7,])))
print(paste0("how many edges (> 0.8): ", nrow(dat[dat$prob>=0.8,])))
print(paste0("how many nodes: ", length(unique(c(as.vector(dat$pre), as.vector(dat$nex))))))
print(paste0("how many cell types: ", length(unique(c(as.vector(dat$pre_cell), as.vector(dat$nex_cell))))))

# extract edges with prob > edge_prob_thresh ===================================

x = dat[dat$prob>=edge_prob_thresh,]
x = x[,c("pre","nex","prob")]
print(paste0("how many nodes now: ", length(unique(c(as.vector(x$pre), as.vector(x$nex))))))
res <- x


dat_sub = res
dat_sub$pre_cell = unlist(lapply(as.vector(dat_sub$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat_sub$nex_cell = unlist(lapply(as.vector(dat_sub$nex), function(x) strsplit(x,"[:]")[[1]][2]))

# summary on edges used to create the tree =====================================

print(paste0("how many edges: ", nrow(res)))
print(paste0("how many nodes: ", length(unique(c(as.vector(res$pre), as.vector(res$nex))))))

nex_list = as.vector(unique(res$nex))
tree = NULL
for(i in 1:length(nex_list)){
  res_sub = res[res$nex==nex_list[i],]
  # if(nrow(res_sub)==1){
    tree = rbind(tree, res_sub)
  # } else {
  #   res_sub = res_sub[order(res_sub$prob, decreasing = TRUE),]
  #   tree = rbind(tree, res_sub[1,])
  # }
}
tree = data.frame(tree)
tree = tree[,c("pre","nex")]

write.table(res, paste0("edge_prob.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(tree, paste0("tree_edge.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
