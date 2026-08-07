#! /usr/bin/env nextflow

log.info """
	TOME timecourse trajectory pipeline
	===================================
	samplesheet: ${params.samplesheet}
	results_dir: ${params.results_dir}
	"""
	.stripIndent()

/* input files:	
 * samplesheet
 * integrated seurat object
 */
 


/*
 * prepare and validate input seurat object
 */

process validate_seurat {
	label 'process_medium'
	container = "tjmgibson/scrnaseq_preprocess:v2"
	
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
	container = "tjmgibson/scrnaseq_preprocess:v2"
	publishDir "${params.results_dir}/pairwise_integration/", mode: 'copy'
	
	input:
	tuple val(timepoint_pair), val(time_1), val(time_2)
	path(seurat_obj)
	path(lineage_knn_source)
	
	output:
	tuple val(timepoint_pair), val(time_1), val(time_2), path("${time_1}_${time_2}.rds"), path("${time_1}_${time_2}_umap3.rds"), path("${time_1}_${time_2}_anno.rds"), path("${time_1}_${time_2}_Knn_umap.rds")
	
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

//  process permute_lineage

//  process permute_lineage {
// 	tag "$timepoint_pair"
// 	label 'process_medium'
// 	container = "tjmgibson/scrnaseq_preprocess:v2"
// 	publishDir "${params.results_dir}/lineage_permutation/", mode: 'copy'
	
// 	input:
// 	tuple val(timepoint_pair), val(time_1), val(time_2), path(integrated_pair_seurat), path(embeddings), path(annotation), path(knn_umap)
	
// 	output:
// 	path("${time_1}_${time_2}_Knn_umap_permutation.rds")
	
// 	script:
//     """
// 	permute_Knn.R ${time_1} ${time_2} ${embeddings} ${annotation} ${params.n_permutations} ${params.k_neighbors}
// 	"""
    
//     stub:
//     """
//     touch ${time_1}_${time_2}_Knn_umap_permutation.rds
//     """	
// }

//  process summarize_results
 process summarize_results {
	label 'process_medium'
	container = "tjmgibson/scrnaseq_preprocess:v2"
	publishDir "${params.results_dir}/tree_results/", mode: 'copy'
	
	input:
	path(knn_umap_files)
	
	output:
	tuple path("edge_all.rds"), path("tree_edge.txt"), path("edge_prob.txt")
	
	script:
    """
	summarize_results.R  ${params.timepoint_names} ${params.n_repeats} ${params.min_edge_weight}
	"""
    
    stub:
    """
    touch edge_all.rds
	touch edge_prob.txt
	touch tree_edge.txt
    """	
}

// // generate json file for tree 
// process create_map {
// 	label 'process_medium'
// 	container = "tjmgibson/scrnaseq_preprocess:v2"
// 	publishDir "${params.results_dir}/tree_results/", mode: 'copy'

// 	input:
// 	tuple path(edge_all), path(tree_edge), path(edge_prob)
// 	path(celltype_group)
	
// 	output:
// 	path("tree.json")
	
// 	script:
//     """
// 	create_map.py  ${params.timepoint_names} $tree_edge $edge_prob ${celltype_group}
// 	"""
    
//     stub:
//     """
//     touch tree.json
//     """	
// }

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

	// permutation_ch = permute_lineage(lineage_ch)
	// | collect

	summary_ch = summarize_results(permutation_ch)

	map_ch = create_map(
		summary_ch,
		file(params.celltype_group, checkIfExists: true)
	)

}