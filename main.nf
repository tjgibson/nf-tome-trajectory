#! /usr/bin/env nextflow


 
/*
 * prepare and validate input seurat object
 */

process validate_seurat {
	label 'process_medium'
	container "tjmgibson/scrnaseq_preprocess:v2"

	input:
	path(seurat_obj)

	output:
	path(seurat_obj)

	script:
	"""
	validate_seurat_obj.R ${seurat_obj} ${params.timepoint_metadata_key} ${params.timepoint_names}
	"""
}

/*
 * compute the lineage
 */

 process compute_lineage {
	tag "$timepoint_pair"
	label 'process_medium'
	container "tjmgibson/scrnaseq_preprocess:v2"
	publishDir "${params.results_dir}/pairwise_integration/", mode: 'copy'
	
	input:
	tuple val(timepoint_pair), val(time_1), val(time_2)
	path(seurat_obj)
	path(lineage_knn_source)
	
	output:
	tuple path("${time_1}_${time_2}.rds"), path("${time_1}_${time_2}_umap3.rds"), path("${time_1}_${time_2}_anno.rds")
	path("${time_1}_${time_2}_Knn_umap.rds"), emit: knn_umap

	script:
	"""
	lineage_connection.R ${seurat_obj} ${params.timepoint_metadata_key} ${params.cell_type_key} ${time_1} ${time_2} ${lineage_knn_source} ${params.k_neighbors}
	"""
    
	stub:
	"""
	touch ${time_1}_${time_2}.rds
	touch ${time_1}_${time_2}_umap3.rds
	touch ${time_1}_${time_2}_anno.rds
	touch ${time_1}_${time_2}_Knn_umap.rds
	"""	
}


//  process summarize_results
 process summarize_results {
	label 'process_low'
	container "tjmgibson/scrnaseq_preprocess:v2"
	publishDir "${params.results_dir}/tree_results/", mode: 'copy'
	
	input:
	path(knn_umap_files)
	
	output:
	tuple path("tree_edge.txt"), path("edge_prob.txt"), path("edge_all.rds")
	path("edge_all.rds"), emit: edge_all
	
	script:
    """
	summarize_results.R  ${params.timepoint_names} ${params.n_repeats} ${params.min_edge_weight}
	"""
    
    stub:
    """
    touch edge_all.rds
	touch edge_prob.txt
	touch tree_edge.txt
	touch tree.pdf
    """	
}

process plot_tree {
	label 'process_low'
	container "tjmgibson/r_networks:v1"
	publishDir "${params.results_dir}/tree_results/", mode: 'copy'

	input:
	path(edge_all)

	output:
	path("tree.pdf")

	script:
    """
	plot_trajectory_tree.R  ${params.timepoint_names} ${params.min_edge_weight} ${edge_all}
	"""
    
    stub:
    """
    touch tree.pdf
    """	

}

/*
 * Run workflow
 */
 
workflow {
	
	validated_seurat_ch = validate_seurat(
		file("${params.merged_seurat_object}",  checkIfExists: true)
		)

	// create channel with adjacent timepoint pairs
	timepoint_list = params.timepoint_names.tokenize(",")
	n_tp_i = timepoint_list.size() - 1

	tp_1_list = timepoint_list[0..n_tp_i-1]
	tp_2_list = timepoint_list[1..n_tp_i]

	tp_pair_list = [tp_1_list, tp_2_list].transpose()
	
	
	tp_pair_ch = Channel.fromList(tp_pair_list)
	| map { time_1, time_2->
        def key = "${time_1}_${time_2}"
        tuple(key, time_1, time_2)
		}


	lineage_ch = compute_lineage(
		tp_pair_ch,
		validated_seurat_ch,
		file(params.lineage_knn_source)
	)
	
	summary_ch = lineage_ch.knn_umap
	| collect
	| summarize_results



	tree_fig_ch = summary_ch.edge_all 
	| plot_tree



}