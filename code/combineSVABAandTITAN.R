#' combineSVABAandTITAN.R
#' author: Gavin Ha 
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date:	  May 7, 2019
#' description: Overlap CNA boundaries to rescue filtered events. Classify SV events.


library(optparse)
option_list <- list(
	make_option(c("--id"), type = "character", help = "Sample ID"),
	make_option(c("--svaba_funcs"), type = "character", help = "Path to file containing SVABA R functions to source."),
	make_option(c("--svabaVCF"), type="character", help = "Path to SVABA barcode rescue VCF file."),
	make_option(c("--titanBinFile"), type="character", help = "Path to TITAN titan.txt output file."),
	make_option(c("--titanSegFile"), type="character", help = "Path to TITAN segs.txt output file."),
	make_option(c("--genomeBuild"), type="character", default="hg19", help = "Genome build: hg19 or hg38. Default [%default]"),
	make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to analyze; string [Default: %default"),
	make_option(c("--minSPAN"), type = "numeric", default = 10000, help = "Minimum span (base pairs) to consider for SV events; does not apply to inversions."),
	make_option(c("--minInversionSPAN"), type = "numeric", default = 1000, help = "Minimum span (base pairs) to consider for inversion events."),
	make_option(c("--outDir"), type="character", help="Path to output directory."),
	make_option(c("--outputSVFile"), type="character", help="Path to output SV file with new annotations."),
	make_option(c("--outputBedpeFile"), type="character", help="Path to output SV file with new annotations as BEDPE."),
	make_option(c("--outputCNFile"), type="character", help="Path to output CNA file with new annotations.")
)

library(data.table)
library(GenomicRanges)
library(VariantAnnotation)
library(stringr)
library(ggplot2)
library(plyr)

options(stringsAsFactors=F, width=165, scipen=999)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

# source utils.R files
source(opt$svaba_funcs)

tumId <- opt$id
svabaVCF <- opt$svabaVCF
cnFile <- opt$titanBinFile
segFile <- opt$titanSegFile
outputSVFile <- opt$outputSVFile
outputCNFile <- opt$outputCNFile
outputBedpeFile <- opt$outputBedpeFile
outDir <- opt$outDir
dir.create(outDir)
outImage <- paste0(outDir, "/", tumId, ".RData")
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
chrs <- as.character(eval(parse(text = opt$chrs)))
minSPAN <- opt$minSPAN
minInvSPAN <- opt$minInversionSPAN

save.image(outImage)

bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
	seqinfo <- Seqinfo(genome=genomeBuild)
} else {
	seqinfo <- seqinfo(get(bsg))
}
seqlevelsStyle(chrs) <- genomeStyle


loess.span <- 0.3
minTrans <- Inf #1e7 #SV size greater than this are considered translocations
numCores <- 6
minBSDsupport <- 4
se.level <- 0.95
filter.quantile <- 0.95
pValThreshold = 1e-10
filterFlags <- c("PASS")
minOverlapFrac <- 0.5
cn.buffer <- 50000 ## for retrieving sv and cn boundary overlaps
sv.cn.change.buffer <- 25000 # for determining prev and next cn segment status at SV
sv.buffer <- 5000 # for determining overlap of sv between tools
maxDupLength <- 10e7 # for determining the max length of tandem duplication
dupSV.bpDiff <- 1000 # for determining duplicate breakpoints in svaba
seg.buffer <- 150000 ## (annotateSVbetweenBkptsWithCN) for determine seg boundary inside of SV breakpoints and minimum seg length to use to look at cn between breakpoints in an event
cn.buffer <- 1e6 ## (getSegSVoverlap) for retrieving sv and cn boundary overlaps for del, gains, inv
cn.invDup.buffer <- cn.buffer * 10 ## (getSegSVoverlap) for sv and cn boundaries overlaps for inverted dups
sv.cn.change.buffer <- 25000 # for determining prev and next cn segment status at SV
minOverlapFrac <- 0.75 # frac overlap of SV with CN during boundary overlap (getSegSVoverlap)
maxDupLength <- 10e7
minNumSeg.simpleSV <- 5
maxInvSPAN <- 5e6
maxFBISPAN <- 30000 # used to determine BX counts - 99% quantile used
minFBIcn <- 4

save.image(outImage)

## load sample list
# rownames = tumor id, column.1=V5 = normal id
#samples <- read.delim(sampleFile, header=F, row.names=1, as.is=T)[,c(4),drop=F]
#colnames(samples) <- c("Tumor", "Normal")

######################################
############# LOAD SVABA #############
######################################
message("Loading svaba results: ", svabaVCF)
svaba <- loadVCFtoDataTableByChromosome(svabaVCF, chr=chrs, genomeStyle=genomeStyle, minSPAN=minSPAN, minBSDsupport=minBSDsupport, dupSV.bpDiff=dupSV.bpDiff)
svaba <- cbind(Sample = tumId, SV.id = 1:nrow(svaba), svaba)

## get SV events that pass filter - intra and interchr events ##
# fold-back inversion 
maxfoldBackInvAD <- svaba[SPAN > 0 & SPAN < maxFBISPAN & orient_1==orient_2, quantile(AD, 0.95)]
indFBI <- svaba[support=="SVABA" & 
			(SPAN > 0 & SPAN < maxFBISPAN & orient_1 == orient_2), #& #fold-back inv
			which=TRUE]
# use minSPAN or maxFBISPAN
indSVABA <- svaba[support=="SVABA" & ((SPAN >= minSPAN ) | SPAN == -1), which=TRUE]
save.image(outImage)

###########################################
###  COPY NUMBER RESCUE OF SVABA EVENTS ###
###########################################
message("Copy number rescue for svaba results...")
segs <- fread(segFile)
segs <- cbind(SEG.id = 1:nrow(segs), segs)	
sv.seg <- getSegSVoverlap(segs, svaba, event.cn=unique(segs$Corrected_Call), buffer=cn.buffer)
sv.seg.interChr <- getSegSVoverlap(segs, svaba, event.cn=unique(segs$Corrected_Call), buffer=cn.buffer, interChr=TRUE)
sv.seg <- rbind(sv.seg, sv.seg.interChr)
# intra-chr
indCN <- sv.seg[support %in% c("SVABA"), SV.id]
svaba[SV.id %in% indCN, support := "CN"]
svaba[!SV.id %in% c(indSVABA, indFBI, indCN), support := NA]
svabaAll <- copy(svaba); setkey(svabaAll, SV.id)
svaba <- svaba[!is.na(support)]
save.image(outImage)


## combine to save data.table ##
svaba[, type := getSVType(svaba, minColSPAN = minSPAN, minTrans = minTrans)]#, maxInvSPAN = maxInvSPAN, maxFBISPAN = maxFBISPAN)]

## output to files ##
outFile <- paste0(outDir, "/", tumId, "_svaba.txt")
write.table(svaba, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

######################################
########### ANNOTATE WITH CN ##########
######################################
message("Annotating SVs with copy number...")
## annotate segment CN overlapping exactly at breakpoints
sv <- copy(svaba)
annot <- annotateSVwithCN(sv, segs, cnColToAnnotate = "Corrected_Copy_Number")
annotMaj <- annotateSVwithCN(sv, segs, cnColToAnnotate = "MajorCN")
annotMin <- annotateSVwithCN(sv, segs, cnColToAnnotate = "MinorCN")
sv[annot$ind1, Copy_Number_1 := annot$annot1]
sv[annot$ind2, Copy_Number_2 := annot$annot2]
sv[annotMaj$ind1, MajorCN_1 := annotMaj$annot1]
sv[annotMin$ind1, MinorCN_1 := annotMin$annot1]
sv[annotMaj$ind2, MajorCN_2 := annotMaj$annot2]
sv[annotMin$ind2, MinorCN_2 := annotMin$annot2]

cn <- fread(cnFile)
cn <- unique(cn[, .(Sample, Chr, Start, End, LogRatio, logR_Copy_Number, Corrected_Copy_Number, Corrected_Call)])
cn <- cn[!is.na(Start) & !is.na(End)]
cn <- cbind(CN.id = 1:nrow(cn), cn)	
sv[, Ploidy.medianCN := median(cn$Corrected_Copy_Number)]
## annotate bin CN to prev/next of both breakpoints
cn[, Corrected_Copy_Number := round(logR_Copy_Number)]
cn[, Corrected_Copy_Number_prev := c(NA, Corrected_Copy_Number[-.N]), by = Chr]
cn[, Corrected_Copy_Number_next := c(Corrected_Copy_Number[-1], NA), by = Chr]
annot <- annotateSVwithCN(sv, cn, cnColToAnnotate = "Corrected_Copy_Number")
sv[annot$ind1, Copy_Number_1 := as.integer(annot$annot1)]
sv[annot$ind2, Copy_Number_2 := as.integer(annot$annot2)]
annot <- annotateSVwithFlankCN(sv, cn, cnColToAnnotate = "Corrected_Copy_Number")
sv[, Copy_Number_1_prev := annot$annot$Corrected_Copy_Number_1_prev]
sv[, Copy_Number_1_next := annot$annot$Corrected_Copy_Number_1_next]
sv[, Copy_Number_2_prev := annot$annot$Corrected_Copy_Number_2_prev]
sv[, Copy_Number_2_next := annot$annot$Corrected_Copy_Number_2_next]
## annotate mean CN across bins between intra-chr SV breakpoints; also # of segments 
annot <- annotateSVbetweenBkptsWithCN(sv, cn, segs, buffer = seg.buffer, 
	cnColToAnnotate = "Corrected_Copy_Number", fun = "mean")
sv[annot$ind, Copy_Number_1_2_mean := round(annot$annot.cn$cn)]
sv[annot$ind, Copy_Number_1_2_numSegs := annot$annot.seg$NumSeg]

sv[, type := getSVType(sv, minColSPAN = minSPAN, minTrans = minTrans)]
save.image(outImage)

########################################################################
######################### SV CLASSIFICATIONS ###########################
########################################################################
if (segs[!grepl("X", Chromosome), median(Corrected_Copy_Number, na.rm=T)] == 2 || 
		segs[!grepl("X", Chromosome), median(Corrected_Copy_Number, na.rm=T)] == 3){
	del.cn <- c("HOMD", "DLOH", "NLOH", "ALOH")
	neut.cn <- c("NEUT", "HET")
}else if (segs[!grepl("X", Chromosome), median(Corrected_Copy_Number, na.rm=T)] >= 4){
	del.cn <- c("HOMD", "DLOH", "NLOH", "ALOH", "GAIN")
	neut.cn <- c("NEUT", "HET", "BCNA")
}
gain.cn <- c("GAIN", "NLOH", "ALOH", "ASCNA", "BCNA", "UBCNA", "AMP", "HLAMP")
#neut.cn <- c("NEUT", "HET")

segs[, TDflankInd := getTDcnFlank(segs, maxDupLength = maxDupLength)]

save.image(outImage)

## events seen in 2 or more tools
#consensus.sv.id <- sv[which(rowSums(!is.na(sv[, .(overlap.SVABA.id, overlap.GROCSVS.id, overlap.LONGRANGER.id)])) >=2), SV.id]
consensus.sv.id <- sv[, SV.id]

########################################################################
## interchromosomal translocations - balanced vs unbalanced ##
########################################################################
#interchr.seg.sv <- getSegSVoverlap(segs, sv, event.cn=neut.cn, buffer=cn.buffer, interChr = TRUE)
interchrUBal.cn.sv <- sv[type == "InterChr" & 
		((Copy_Number_1_prev != Copy_Number_1_next) | (Copy_Number_2_prev != Copy_Number_2_next))]
		#((orient_1 == "rev" & Copy_Number_1_prev > Copy_Number_1_next) | 
		# (orient_1 == "fwd" & Copy_Number_1_prev < Copy_Number_1_next)) &
		#((orient_2 == "rev" & Copy_Number_2_prev > Copy_Number_2_next) |
		# (orient_2 == "fwd" & Copy_Number_2_prev < Copy_Number_2_next)) |
if (nrow(interchrUBal.cn.sv) > 0){
	interchrUBal.sv.id <- sort(interchrUBal.cn.sv$SV.id)
	sv[SV.id %in% interchrUBal.sv.id, CN_overlap_type := "Trans-Unbal"]
	interchrBAL.cn.sv <- sv[type == "InterChr" &
			(Copy_Number_1_prev == Copy_Number_1_next) & (Copy_Number_2_prev == Copy_Number_2_next)]
	interchrBAL.sv.id <- sort(interchrBAL.cn.sv$SV.id)
	sv[SV.id %in% interchrBAL.sv.id, CN_overlap_type := "Trans-Bal"]
	# exclude BX rescue for interchr if not seen by another tool
	sv[type == "InterChr" & !SV.id %in% consensus.sv.id & (support == "BX"), CN_overlap_type := NA] 
}

########################################################################
## simple deletions ##
########################################################################
# use segment overlap and
# bin-level adjacent CN: b1.prev > b1.next & b2.prev < b2.next & meanCN between b1-b2 < b1.prev/b2.next 
del.cn.sv <- sv[type == "Deletion" & 
	((Copy_Number_1_prev > Copy_Number_1_next | Copy_Number_2_prev < Copy_Number_2_next) &
	(Copy_Number_1_prev > Copy_Number_1_2_mean | Copy_Number_2_next > Copy_Number_1_2_mean)) &
	(Copy_Number_1_2_numSegs <= minNumSeg.simpleSV)]
del.seg.sv <- getSegSVoverlap(segs, sv, event.cn=del.cn, buffer=cn.buffer)	
if (!is.null(del.seg.sv)){
del.sv <- del.seg.sv[type == "Deletion" &
	sv.overlap.frac >= minOverlapFrac & 
	bkpt1.1 < cn.buffer*5 & bkpt2.2 < cn.buffer*5 &
	(Copy_Number_1_prev > Copy_Number_1_2_mean & Copy_Number_2_next > Copy_Number_1_2_mean)]
}else{ del.sv <- NULL }
del.sv.id <- sort(union(del.sv$SV.id, del.cn.sv$SV.id))
sv[SV.id %in% del.sv.id, CN_overlap_type := "Deletion"]
segs[SEG.id %in% del.sv$SEG.id, SV_overlap_type := "Deletion"]
if (!is.null(del.sv)){
	del.sv.seg.id <- del.sv[, paste0(SEG.id, collapse=";"), by=SV.id][order(SV.id)]
	sv[SV.id %in% del.sv.seg.id$SV.id, SEG.overlap.id := del.sv.seg.id$V1]
	del.seg.sv.id <- del.sv[, paste0(SV.id, collapse=";"), by=SEG.id][order(SEG.id)]
	segs[SEG.id %in% del.seg.sv.id$SEG.id, SV.overlap.id := del.seg.sv.id$V1]
}
sv[is.na(CN_overlap_type) & type=="Deletion" & SPAN > minSPAN & SPAN < cn.buffer & Copy_Number_1_2_numSegs <= 2, CN_overlap_type := "Deletion"]

########################################################################
## simple tandem duplications ##
########################################################################
# use segment overlap and breakpoint CN
gain.cn.sv <- sv[type == "Duplication" & 
	((Copy_Number_1_prev < Copy_Number_1_next | Copy_Number_2_prev > Copy_Number_2_next) &
	(Copy_Number_1_prev < Copy_Number_1_2_mean | Copy_Number_2_next < Copy_Number_1_2_mean)) &
	(Copy_Number_1_2_numSegs <= minNumSeg.simpleSV)]
# bin-level adjacent CN: b1.prev < b2.next & b2.prev > b2.next & meanCN between b1-b2 > b1.prev/b2.next
gain.seg.sv <- getSegSVoverlap(segs[TDflankInd==TRUE], sv, event.cn=gain.cn, buffer=cn.buffer)	
if (!is.null(gain.seg.sv)){
	gain.sv <- gain.seg.sv[type == "Duplication" & 
		sv.overlap.frac >= minOverlapFrac & 
		bkpt1.1 < cn.buffer*5 & bkpt2.2 < cn.buffer*5 &
		(Copy_Number_1_prev < Copy_Number_1_2_mean | Copy_Number_2_next < Copy_Number_1_2_mean)]
}else{ gain.sv <- NULL }
gain.sv.id <- sort(union(gain.sv$SV.id, gain.cn.sv$SV.id))	
sv[SV.id %in% gain.sv.id, CN_overlap_type := "TandemDup"]
segs[gain.sv$SEG.id, SV_overlap_type := "TandemDup"]
segs[gain.sv$SEG.id, SV_start_1 := gain.sv$start_1]
segs[gain.sv$SEG.id, SV_start_2 := gain.sv$start_2]
#sv[SV.id %in% gain.sv.id, Copy_Number_SVoverlap := gain.sv$Corrected_Copy_Number]
#sv[SV.id %in% gain.sv.id, Copy_Number_Call := gain.sv$Corrected_Call]
#sv[SV.id %in% gain.sv.id, Median_HaplotypeRatio := gain.sv$Median_HaplotypeRatio]
if (!is.null(gain.sv)){
	gain.sv.seg.id <- gain.sv[, paste0(SEG.id, collapse=";"), by=SV.id][order(SV.id)]
	sv[SV.id %in% gain.sv.seg.id$SV.id, SEG.overlap.id := gain.sv.seg.id$V1]
	gain.seg.sv.id <- gain.sv[, paste0(SV.id, collapse=";"), by=SEG.id][order(SEG.id)]
	segs[SEG.id %in% gain.seg.sv.id$SEG.id, SV.overlap.id := gain.seg.sv.id$V1]
}
#sv[is.na(CN_overlap_type) & type=="Duplication" & Tool != "LONGRANGER" & SPAN > minSPAN & SPAN < cn.buffer & Copy_Number_1_2_numSegs <= 2, CN_overlap_type := "TandemDup"]
sv[is.na(CN_overlap_type) & type=="Duplication" & SPAN > minSPAN & SPAN < cn.buffer & Copy_Number_1_2_numSegs <= 2, CN_overlap_type := "TandemDup"]

########################################################################
## unbalanced inversions ##
########################################################################
inv.cn.sv <- sv[type == "Inversion" & SPAN >= minInvSPAN & SPAN <= maxInvSPAN &
	((Copy_Number_1_prev != Copy_Number_1_2_mean | Copy_Number_2_next != Copy_Number_1_2_mean) |
	(Copy_Number_1_prev != Copy_Number_1_next | Copy_Number_2_prev != Copy_Number_2_next))]
#inv.seg.sv <- getSegSVoverlap(segs, sv, event.cn=neut.cn, buffer=cn.buffer)
#if (!is.null(inv.seg.sv)){
#inv.sv <- inv.seg.sv[type == "Inversion" & 
#	sv.overlap.frac >= minOverlapFrac & 
#	bkpt1.1 < cn.buffer & bkpt2.2 < cn.buffer &
#	(Copy_Number_1_prev != Copy_Number_1_2_mean & Copy_Number_2_next != Copy_Number_1_2_mean)]
#}else{ inv.sv <- NULL }
#inv.sv.id <- sort(union(inv.sv$SV.id, inv.cn.sv$SV.id))	
sv[SV.id %in% inv.cn.sv$SV.id, CN_overlap_type := "Inversion-Unbal"]
#segs[SEG.id %in% inv.sv$SEG.id, SV_overlap_type := "Inversion-UnBal"]
#if (!is.null(inv.sv)){
#	inv.sv.seg.id <- inv.sv[, paste0(SEG.id, collapse=";"), by=SV.id][order(SV.id)]
#	sv[SV.id %in% inv.sv.seg.id$SV.id, SEG.overlap.id := inv.sv.seg.id$V1]
#	inv.seg.sv.id <- inv.sv[, paste0(SV.id, collapse=";"), by=SEG.id][order(SEG.id)]
#	segs[SEG.id %in% inv.seg.sv.id$SEG.id, SV.overlap.id := inv.seg.sv.id$V1]
#}

########################################################################
## balanced inversions ##
########################################################################
inv.cn.sv <- sv[type == "Inversion" & SPAN >= minInvSPAN & SPAN <= maxInvSPAN &
	((Copy_Number_1_prev == Copy_Number_1_2_mean & Copy_Number_2_next == Copy_Number_1_2_mean) |
	(Copy_Number_1_prev == Copy_Number_1_next & Copy_Number_2_prev == Copy_Number_2_next))]
#inv.seg.sv <- getSegSVoverlap(segs, sv, event.cn=neut.cn, buffer=cn.buffer)
#if (!is.null(inv.seg.sv)){
#inv.sv <- inv.seg.sv[type == "Inversion" & 
#	sv.overlap.frac >= minOverlapFrac & 
#	bkpt1.1 < cn.buffer & bkpt2.2 < cn.buffer &
#	(Copy_Number_1_prev == Copy_Number_1_2_mean & Copy_Number_2_next == Copy_Number_1_2_mean)]
#}else{ inv.sv <- NULL }
#inv.sv.id <- sort(union(inv.sv$SV.id, inv.cn.sv$SV.id))	
sv[SV.id %in% inv.cn.sv$SV.id, CN_overlap_type := "Inversion-Bal"]
#segs[SEG.id %in% inv.sv$SEG.id, SV_overlap_type := "Inversion-Bal"]
#if (!is.null(inv.sv)){
#	inv.sv.seg.id <- inv.sv[, paste0(SEG.id, collapse=";"), by=SV.id][order(SV.id)]
#	sv[SV.id %in% inv.sv.seg.id$SV.id, SEG.overlap.id := inv.sv.seg.id$V1]
#	inv.seg.sv.id <- inv.sv[, paste0(SV.id, collapse=";"), by=SEG.id][order(SEG.id)]
#	segs[SEG.id %in% inv.seg.sv.id$SEG.id, SV.overlap.id := inv.seg.sv.id$V1]
#}

########################################################################
## foldback inversions ##
########################################################################
inv.cn.sv <- sv[type == "Inversion" & (SPAN < maxFBISPAN) &
	(Copy_Number_1_prev != Copy_Number_2_next | Copy_Number_1 != Copy_Number_2) &
	(Copy_Number_1 / Ploidy.medianCN * 2 > minFBIcn | Copy_Number_2 / Ploidy.medianCN * 2 > minFBIcn)]
sv[SV.id %in% inv.cn.sv$SV.id, CN_overlap_type := "Inversion-FoldBack"]

#sv[type == "FoldBackInv-H", CN_overlap_type := "FoldBackInv-H"]
#sv[type == "FoldBackInv-T", CN_overlap_type := "FoldBackInv-T"]

########################################################################
## Rest: balanced rearrangements - intra-chr SVs ##
########################################################################
sv[is.na(CN_overlap_type) & type=="Inversion" & SPAN >= minInvSPAN &
	((Copy_Number_1_prev == Copy_Number_1_2_mean & Copy_Number_2_next == Copy_Number_1_2_mean) |
	(Copy_Number_1_prev == Copy_Number_1_next & Copy_Number_2_prev == Copy_Number_2_next)),
	CN_overlap_type := "Balanced"]
## Rest: unbalanced rearrangements ##
sv[is.na(CN_overlap_type) & !is.na(type) & SPAN >= minInvSPAN, CN_overlap_type := "Unbalanced"]

########################################################################
## filter short inversions and deletions (longranger & svaba) ##
########################################################################
#sv <- sv[!is.na(CN_overlap_type) | !is.na(overlap.GROCSVS.id)]
sv <- sv[!is.na(CN_overlap_type)]

########################################################################
## inverted duplications (templated insertions?) ##
########################################################################
# use segment overlap and 
# bin-level adjacent CN: b1.prev < b2.next & b2.prev > b2.next & meanCN between b1-b2 > b1.prev/b2.next
#sv.id.toUse <- sv[is.na(CN_overlap_type), SV.id]
#gain.seg.sv <- getSegSVoverlap(segs, sv[Tool != "LONGRANGER" & SV.id %in% sv.id.toUse], event.cn=gain.cn, buffer=cn.invDup.buffer)
#gain.cn.sv <- sv[Tool != "LONGRANGER" & type == "Inversion" & 
#		(Copy_Number_1_prev < Copy_Number_1_next & Copy_Number_2_prev > Copy_Number_2_next) &
#		(Copy_Number_1_prev < Copy_Number_1_2_mean & Copy_Number_2_next < Copy_Number_1_2_mean) &
#		(abs(Copy_Number_1_next - Copy_Number_2_prev) == 0 | Copy_Number_1_2_numSegs > 2)]
#bkpt.a.b - a is sv breakpoint 1/2 and b is seg start/end
#if (!is.null(gain.seg.sv)){
#gain.sv <- gain.seg.sv[#SPAN > maxInvSPAN &
#	xor(xor(bkpt1.1 < cn.buffer & Copy_Number_1_prev < Copy_Number_1_next & orient_1 == "fwd",
#			 bkpt2.1 < cn.buffer & Copy_Number_2_prev < Copy_Number_2_next & orient_2 == "fwd"),
#	 xor(bkpt1.2 < cn.buffer & Copy_Number_1_prev > Copy_Number_1_next & orient_1 == "rev",
#			 bkpt2.2 < cn.buffer & Copy_Number_2_prev > Copy_Number_2_next & orient_2 == "rev"))]
#}else{ gain.sv <- NULL}
#tempIns.seg.id <- sort(unique(gain.sv[, .N > 1, by=SEG.id][V1==TRUE, SEG.id]))
#tempIns.sv.id <- sort(unique(gain.sv[SEG.id %in% tempIns.seg.id, SV.id])
#sv[SV.id %in% tempIns.sv.id, CN_overlap_type := "TempIns"]
#segs[SEG.id %in% tempIns.seg.id, SV_overlap_type := "TempIns"]	
#if (!is.null(gain.sv)){
#	gain.sv.seg.id <- gain.sv[, paste0(SEG.id, collapse=";"), by=SV.id][order(SV.id)]
#	sv[SV.id %in% gain.sv.seg.id$SV.id, SEG.overlap.id := gain.sv.seg.id$V1]
#	gain.seg.sv.id <- gain.sv[, paste0(SV.id, collapse=";"), by=SEG.id][order(SEG.id)]
#	segs[SEG.id %in% gain.seg.sv.id$SEG.id, SV.overlap.id := gain.seg.sv.id$V1]
#}

########################################################################
## output files ##
########################################################################
save.image(outImage)
#sort by chr
sv$chromosome_1 = factor(sv$chromosome_1, levels=chrs)
sv$chromosome_2 = factor(sv$chromosome_2, levels=chrs)
sv <- sv[order(chromosome_1, chromosome_2, start_1, start_2)]
write.table(sv, file=outputSVFile, col.names=T, row.names=F, quote=F, sep="\t")
write.table(segs, file=outputCNFile, col.names=T, row.names=F, quote=F, sep="\t")


########################################################################
## output bedpe files ##
########################################################################
sv.out <- copy(sv)
sv.out[, c("chrom1", "start1", "end1", "chrom2", "start2", "end2") := .(chromosome_1, start_1, start_1, chromosome_2, start_2, start_2)]
sv.out[, strand1 := { strand = rep("+", .N); strand[orient_1=="fwd"] = "-"; strand }]
sv.out[, strand2 := { strand = rep("+", .N); strand[orient_2=="fwd"] = "-"; strand }]
sv.out[, name := "."]; sv.out[, score := "."]
colNameOrder <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2")
setcolorder(sv.out, c(colNameOrder, colnames(sv.out)[!colnames(sv.out) %in% colNameOrder]))
write.table(sv.out, file = outputBedpeFile, col.names=T, row.names=F, quote=F, sep="\t", na=".")


