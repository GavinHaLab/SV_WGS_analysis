#' create_circosPlot.R
#' author: Minjeong Ko
#' Fred Hutchinson Cancer Research Center
#' contact: <mko@fredhutch.org>
#' date:  March 19, 2020
#' description: Generate circos plots of copy number and SV breakpoint arcs.

library(RCircos)
library(optparse)

option_list <- list(
  make_option(c("--id"), type = "character", help = "Sample ID"),
  make_option(c("--svFile"), type = "character", help = "Combined Svaba and Titan SV file (*svabaTitan.sv.bedpe)."),
  make_option(c("--cnFile"), type = "character", help = "Combined Svaba and Titan CNA file (*svabaTitan.cn.txt)."),  
  make_option(c("--genomeBuild"), type = "character", help = "genomeBuild"),
  make_option(c("--outPlotFile"), type="character", help="Path to output figure file.")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

id <- opt$id
svFile <- opt$svFile
cnFile <- opt$cnFile
genomeBuild <- opt$genomeBuild
outPlotFile <- opt$outPlotFile

#Load chromosome cytoband data
if(genomeBuild =="hg38"){
	data(UCSC.HG38.Human.CytoBandIdeogram)
	cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
	}else{
	data(UCSC.HG19.Human.CytoBandIdeogram)
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
}

#Setup RCircos core components
RCircos.Set.Core.Components(cyto.info, tracks.inside=10, tracks.outside=0) 

sv.data <- read.table(svFile, sep="\t", head=T)
cnv.data <- read.table(cnFile, sep="\t", head=T)[,c("Chromosome","Start","End","Corrected_Copy_Number")]


pdf(outPlotFile,onefile=T,width=19,height=19)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Histogram.Plot(cnv.data, data.col=4, track.num=1, side="in")
RCircos.Link.Plot(sv.data, track.num=2, by.chromosome=FALSE)
dev.off()
