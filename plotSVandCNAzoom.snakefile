configfile: "config/configPlotZoom.yaml"
configfile: "config/samples.yaml"

import glob
def getTITANpath(base, id, ext):
  return glob.glob(''.join([base, "results/titan/hmm/optimalClusterSolution/", id, "_cluster*", ext]))

#TUM, CLUST = glob_wildcards("../../TITAN/snakemake/results/titan/optimalClusterSolution/{tum}_cluster1.titan.ichor.cna.txt")
#SEG,CLUST = glob_wildcards(config["titanPath"], "/results/titan/optimalClusterSolution/{tumor}_cluster{clust}.titan.ichor.cna.txt")


rule all:
  input: 
  	expand("results/plotSVABAandTITAN_zoom/{plotID}/{tumor}_CNA-SV_{type}_chr{chr}-{start}-{end}.{format}", tumor=config["pairings"], plotID=config["plot_id"], type=config["plot_type"], chr=config["plot_chr"], start=config["plot_startPos"], end=config["plot_endPos"], format=config["plot_format"])
	
		
rule plotSVABAandTITAN:
	input:
		svabaVCF="results/combineSVABAandTITAN/{tumor}/{tumor}.svabaTitan.sv.txt",
		#svabaVCF=expand("results/combineSVABAandTITAN/{tumor}/{tumor}.svabaTitan.sv.txt", tumor=config["pairings"]),
		titanBinFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.cna.txt"),
		titanSegFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".titan.ichor.seg.noSNPs.txt"),
		titanParamFile=lambda wildcards: getTITANpath(config["titan_results"], wildcards.tumor, ".params.txt")
	output:
		"results/plotSVABAandTITAN_zoom/{plotID}/{tumor}_CNA-SV_{type}_chr{chr}-{start}-{end}.{format}"
	params:
		plotSVCNscript=config["plotSVCN_script"],
		svabafuncs=config["svaba_funcs"],
		plotfuncs=config["plot_funcs"],
		libdir=config["titan_libdir"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		cytobandFile=config["cytobandFile"],
		zoom=config["plot_zoom"],
		chr=config["plot_chr"],
		start=config["plot_startPos"],
		end=config["plot_endPos"],
		ylim=config["plot_ylim"],
		type=config["plot_type"],
		geneFile=config["plot_geneFile"],
		size=config["plot_size"],
		format=config["plot_format"]	
	log:
		"logs/plotSVABAandTITAN_zoom/{plotID}/{tumor}_{type}_chr{chr}-{start}-{end}.{format}.log"
	shell:
		"Rscript {params.plotSVCNscript} --id {wildcards.tumor} --svaba_funcs {params.svabafuncs} --plot_funcs {params.plotfuncs} --titan_libdir {params.libdir} --svFile {input.svabaVCF} --titanBinFile {input.titanBinFile} --titanSegFile {input.titanSegFile} --titanParamFile {input.titanParamFile} --chrs \"{params.chr}\" --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --cytobandFile {params.cytobandFile} --start {params.start} --end {params.end} --zoom {params.zoom} --plotYlim \"{params.ylim}\" --geneFile {params.geneFile} --plotCNAtype {params.type} --plotSize \"{params.size}\" --outPlotFile {output} > {log} 2> {log}" 
	

