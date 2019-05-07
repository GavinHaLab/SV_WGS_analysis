configfile: "config/config.yaml"
configfile: "config/samples.yaml"


rule svabaAll:
  input: 
  	#expand("results/svaba/{tumor}", tumor=config["pairings"]),
  	expand("results/svaba/{tumor}/{tumor}.svaba.somatic.sv.vcf", tumor=config["pairings"]),
  	#expand("results/svaba/{tumor}/{tumor}.svaba.unfiltered.somatic.sv.vcf", tumor=config["pairings"])

rule runSvaba:
	input:
		tum=lambda wildcards: config["samples"][wildcards.tumor],
		norm=lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
	output:		
		somaticFiltVCF="results/svaba/{tumor}/{tumor}.svaba.somatic.sv.vcf",
		#somaticUnfiltVCF="results/svaba/{tumor}/{tumor}.svaba.unfiltered.somatic.sv.vcf",
		#bps="results/svaba/{tumor}/{tumor}.bps.txt.gz"
	params:
		outDir="results/svaba/{tumor}",
		svabaExe=config["svaba_exe"],
		refGenome=config["refGenome"],
		numThreads=config["svaba_numThreads"],
		dbSNPindelVCF=config["svaba_dbSNPindelVCF"],
		mem=config["svaba_mem"],
		runtime=config["svaba_runtime"],
		pe=config["svaba_numCores"]
	log:
		"logs/svaba/{tumor}.log"
	shell:
		"{params.svabaExe} run -t {input.tum} -n {input.norm} -G {params.refGenome} -p {params.numThreads} -D {params.dbSNPindelVCF} -a {params.outDir}/{wildcards.tumor} > {log} 2> {log}"

	

