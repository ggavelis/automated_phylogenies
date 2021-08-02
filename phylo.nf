#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.unaligned = "/mnt/scgc/scgc_nfs/lab/ggavelis/sandbox/phylo_nx/phylogenetic_pipeline/query/1_to_align.faa"


params.outdir = "results"
params.mafft_options = "--reorder --bl 30 --op 1.0 --maxiterate 1000 --retree 1 --genafpair --quiet"
params.trimal_options = "-gappyout"
params.iqtree_options = "-m TEST -alrt 1000 -bb 1000"
params.iqtree_outgroup = "-o "

Channel
    .fromPath( params.unaligned )
    .set{ fasta_ch }
    
process mafft {
	publishDir params.outdir
	module 'singularity/3.8.0'
	
	input:
	path(unaligned)
	val(params.mafft_options)
	
	output:
	path "2_aligned.faa", emit: aligned
	
	"""
	mafft ${params.mafft_options} ${unaligned}> 2_aligned.faa
	"""
}

process trimal {
	publishDir params.outdir
	module 'singularity/3.8.0'
	
	input:
	path(aligned)
	val(params.trimal_options)
	
	output:
	path "3_trimmed_aligned.faa", emit: trimmed
	
	"""
	trimal -in ${aligned} -out 3_trimmed_aligned.faa ${params.trimal_options}
	"""
}

process iqtree {
	publishDir params.outdir
	module 'singularity/3.8.0'
	cpus 24
	
	input:
	path(trimmed)
	val(params.iqtree_options)
	
	output:
	path "4_iqtree.iqtree", emit: tree_summary
	path "4_iqtree.contree", emit: contree
	path "4_iqtree.log", emit: treelog
	path "4_iqtree.tre", emit: tre
	
	"""
	iqtree -s ${trimmed} ${params.iqtree_options} -pre 4_iqtree -nt ${task.cpus}
	mv 4_iqtree.treefile ./4_iqtree.tre
	"""

}


workflow {
	mafft(fasta_ch, params.mafft_options)
	trimal(mafft.output.aligned, params.trimal_options)
	iqtree(trimal.output.trimmed, params.iqtree_options)
}
