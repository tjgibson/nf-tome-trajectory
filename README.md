---
editor_options: 
  markdown: 
    wrap: sentence
---

## Introduction

This repository contains a nextflow pipeline for performing lineage analysis in single-cell RNA-seq time course data, as described in [Qiu et al., 2022](https://doi.org/10.1038/s41588-022-01018-x).
I adapted the code from their paper ([see github repo here](https://github.com/ChengxiangQiu/tome_code/tree/main)) into a nextflow pipeline in order to make this approach easier to run on new datasets and to make it easier to deploy on a computing cluster.
I only changed code to handle inputs and outputs, provide informative progress messages.
I did not alter any of the code affecting the k-nearest neighbors approach used to infer the ancestors and descendents of a cell type in adjacent time points.

## The original approach

For a detailed description of this method, see the [original paper](https://doi.org/10.1038/s41588-022-01018-x).
Briefly, the authors present a new method called TOME (trajectories of mouse embryogenesis) that can be used for inferring developmental trajectories in single-cell time course data.
The approach performs integration of each pair of adjacent timepoints, then uses a k-nearest neighbors approach to infer which cell type in the adjacent timepoint is closets.
Performing this across a whole timecourse allows connecting celltypes across all time points based on their inferred relationships.

## What this pipeline does

This pipeline uses the minimal amount of code from the original paper to compute the inferred trajectory.
Compared to the original paper, it is worth noting the following points:

-   In the original paper, the authors use both CCA integration and RPCA integration (as implemented in Seurat) to integrate adjacent time points.
    Currently, this pipeline uses only RPCA integration, as it is straightforward and scales well to larger datasets.

-   The original authors use permutation of the connection to compute random background probabilities for each potential trajectory connection.
    As these probabilities are not used for assembling the trajectory tree, I have omitted this step.

-   The original authors use the d3 javascript library for creating a publication-quality developmental trajectory tree.
    I use the networkD3 R library to create a similar tree.
    In the future, I may try to implement using d3.js to create a cleaner figure.

## To-do

-   add options for multiple integration algorithms.
-   Use d3.js for creating a more aesthetic final figure.
