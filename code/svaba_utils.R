###########################################
### Extract SV to data.frame from VCF #####
###########################################
getSVfromCollapsedVCF <- function(vcf, chrs = c(1:22, "X"), genomeStyle = "NCBI"){
	if (class(vcf) == "CollapsedVCF"){
			svGR <- as.data.frame(rowRanges(vcf))
			svGR$ALT <- unlist(svGR$ALT)
			bkpt2 <- str_match(svGR$ALT, "([chr0-9XYA-Z.]+):(\\d+)")[, 2:3]
			mateID <- info(vcf)$MATEID
			sv <- data.frame(chromosome_1=as.vector(svGR$seqnames), start_1=as.numeric(svGR$start), 
																chromosome_2=bkpt2[, 1], start_2=as.numeric(bkpt2[, 2]), REF=svGR$REF,
																mateID=mateID, alt_1=svGR$ALT, alt_2=NA, FILTER=rowRanges(vcf)$FILTER,
																stringsAsFactors=FALSE, orient_1=NA, orient_2=NA)
			## get read support ##
			# BXOL - barcode overlap
			if (sum(rownames(info(header(vcf))) %in% c("BXOL")) > 0){
				sv <- cbind(sv, BXOL = info(vcf)$BXOL)
			}
			# BXC.1 - barcode count for breakpoint 1
			if (sum(rownames(info(header(vcf))) %in% c("BXC.1")) > 0){
				sv <- cbind(sv, BXC.1 = info(vcf)$BXC.1)
			}
			# BXC.2 - barcode count for breakpoint 1
			if (sum(rownames(info(header(vcf))) %in% c("BXC.2")) > 0){
				sv <- cbind(sv, BXC.2 = info(vcf)$BXC.2)
			}
			# BXSupport - barcode count for DR and SR read support
			if (sum(rownames(info(header(vcf))) %in% c("BXSupport")) > 0){
				sv <- cbind(sv, BXSupport = info(vcf)$BXSupport)
			}
			# DR - discordant reads
			if (sum(rownames(geno(header(vcf))) %in% c("DR")) > 0){
				sv <- cbind(sv, DR = geno(vcf)$DR[,ncol(geno(vcf)$DR)])
			}
			# SR - split reads
			if (sum(rownames(geno(header(vcf))) %in% c("SR")) > 0){
				sv <- cbind(sv, SR = geno(vcf)$SR[,ncol(geno(vcf)$SR)])
			}
			# AD - split reads
			if (sum(rownames(geno(header(vcf))) %in% c("AD")) > 0){
				sv <- cbind(sv, AD = geno(vcf)$AD[,ncol(geno(vcf)$AD)])
			}
			# GT - filter out 1/1
			if (sum(rownames(geno(header(vcf))) %in% c("GT")) > 0){
				sv <- cbind(sv, GT = geno(vcf)$GT[,ncol(geno(vcf)$GT)])
			}
			# SPAN 
			if (sum(rownames(info(header(vcf))) %in% c("SPAN")) > 0){
				sv <- cbind(sv, SPAN = info(vcf)$SPAN)
			}
			# SCTG
			if (sum(rownames(info(header(vcf))) %in% c("SCTG")) > 0){
				sv <- cbind(sv, SCTG = info(vcf)$SCTG)
			}
			# BXTags - barcodes for DR and SR read support
			if (sum(rownames(info(header(vcf))) %in% c("BXTags")) > 0){
				sv <- cbind(sv, BXTags = unlist(lapply(as.list(info(vcf)$BXTags), paste, collapse=",")))
				sv[sv$BXTags=="", "BXTags"] <- NA
			}
			## get breakpoint 2 alt ##
			sv$alt_2 <- sv[sv$mateID, "alt_1"]
      
      ## filter by chromosome ##
      		seqlevelsStyle(chrs) <- genomeStyle
      		sv$chromosome_1 <- mapSeqlevels(sv$chromosome_1, style = genomeStyle)
      		sv$chromosome_2 <- mapSeqlevels(sv$chromosome_2, style = genomeStyle)
      		ind <- sv$chromosome_1 %in% chrs & sv$chromosome_2 %in% chrs
			sv <- sv[ind,]
      
			## get orientation ##
			## check the other alt to see how the current breakpoint is oriented
			orient1 <- str_match(sv$alt_2, "\\]|\\[") 
			orient2 <- str_match(sv$alt_1, "\\]|\\[")
			sv$orient_1[orient1 == "]"] <- "rev"
			sv$orient_1[orient1 == "["] <- "fwd"
			sv$orient_2[orient2 == "]"] <- "rev"
			sv$orient_2[orient2 == "["] <- "fwd"
	}else{
		message("vcf object is not a CollapsedVCF")
		sv <- NULL
	}
	# 
	sv <- sv[!is.na(sv$chromosome_1) & !is.na(sv$chromosome_2), ]
	sv <- sv[sv$chromosome_1 %in% chrs & sv$chromosome_2 %in% chrs, ]
	sv <- sv[!is.na(sv$alt_1) & !is.na(sv$alt_2) & !is.na(sv$orient_1) & !is.na(sv$orient_2), ]
	sv <- removeDupSV(sv)
	return(data.table(sv))
}

###############################################
##### set SvABA VCF output genome style #######
###############################################
setVCFgenomeStyle <- function(x, genomeStyle = "NCBI"){
	oldGenomeStyle <- seqlevelsStyle(x)[1]
	if (oldGenomeStyle != genomeStyle){
		message("Changing genomeStyle ", oldGenomeStyle, " to ", genomeStyle, "...")
		# main chromosome field (rowRanges)
		seqlevelsStyle(x) <- genomeStyle
		# need to set the ALT field
		svGR <- as.data.frame(rowRanges(x))
		svGR$ALT <- unlist(svGR$ALT)
		bkpt2 <- str_match(svGR$ALT, "([A-Z\\]\\[]+)([0-9XYA-Zchr.]+):(\\d+)([A-Z\\]\\[]+)")[, -1]
		bkpt2[,2] <- mapSeqlevels(bkpt2[,2], style=genomeStyle)
		svGR$ALT <- paste0(bkpt2[, 1], bkpt2[, 2], ":", bkpt2[, 3], bkpt2[, 4])
		rowRanges(x)$ALT = CharacterList(as.list(svGR$ALT))
	}else{
		message("VCF object is already in ", genomeStyle, " genomeStyle.")
	}
	return(x)
} 

##########################################
##### remove duplicate breakpoints #######
##########################################
removeDupSV <- function(sv){
 #if (nrow(sv) >= 2){
 #		dupInd <- NULL; keepInd <- NULL
#			for (j in 1:nrow(sv)){
#				if (!j %in% dupInd){
#					dupInd <- c(dupInd, which(sv[, "chromosome_2"] %in% sv[j, "chromosome_1"] &
#												sv[, "start_2"] %in% sv[j, "start_1"] &
#												sv[, "start_1"] %in% sv[j, "start_2"]))	
#					keepInd <- c(keepInd, j)
#				}
#			}
	#dupInd <- which(duplicated(str_match(sv$mateID, "\\d+")[,1]))
	#sv <- sv[-dupInd, ]
	sv <- sv[grepl(":1", rownames(sv)), ] ## keep 1st mate in pair 
	## inter chromosomal events #
#	 }
	return(sv)
}

## remove duplicate events if one of them is already PASS #
# can occur if one of start_1 or start_2 is the same
# also remove when coordinates of start_2 are within bp.diff  
removeDupNonPassSV <- function(sv, bp.diff = 5){
  start1DupInd <- sv[duplicated(start_1), , which=T]
  start2DupInd <- sv[duplicated(start_2), , which=T]
  #bothDupInd <- intersect(start1DupInd, start2DupInd)
  bothDupInd <- unique(sort(c(start1DupInd, start2DupInd)))
  removeInd <- NULL
  ## exact duplicates ##
  message("Removing duplicate breakpoints...")
  for (i in 1:length(bothDupInd)){
    chr1 <- sv[bothDupInd[i], chromosome_1]
    start1 <- sv[bothDupInd[i], start_1]
    chr2 <- sv[bothDupInd[i], chromosome_2]
    start2 <- sv[bothDupInd[i], start_2]
    dupInd1 <- sv[chromosome_1 == chr1 & start_1 == start1, , which=T]
    dupInd2 <- sv[chromosome_2 == chr2 & start_2 == start2, , which=T]
    if (length(dupInd1) > 1){ # dup in start_1
      #if there is at least 1 dup that is PASS
      removeInd <- c(removeInd, selectDupSVtoRemove(sv, dupInd1))
    }
    if (length(dupInd2) > 1){ #dup in start_2
      #if there is at least 1 dup that is PASS
      removeInd <- c(removeInd, selectDupSVtoRemove(sv, dupInd2))
    }
  }
  keepInd <- setdiff(1:nrow(sv), removeInd)
  sv <- sv[keepInd]
  ## similar breakpoints within bp.diff for both breakpoints
  message("Removing similar (", bp.diff,"bp) breakpoints ...")
  sortChrInd <- order(sv$chromosome_1, sv$start_1, sv$chromosome_2, sv$start_2)
  diffInd <- sv[sortChrInd, 
    abs(diff(start_1)) <= bp.diff & chromosome_1[-1] == chromosome_1[-nrow(sv)] &
    abs(diff(start_2)) <= bp.diff & chromosome_2[-1] == chromosome_2[-nrow(sv)]]
  diffInd <- sortChrInd[diffInd]
  removeDiffInd <- NULL
  for (i in 1:length(diffInd)){
    diffInd1 <- c(diffInd[i], diffInd[i] + 1)
    if (length(unique((sv[diffInd1, orient_1]))) == 1){ #same orientation
      removeDiffInd <- c(removeDiffInd, selectDupSVtoRemove(sv, diffInd1))
    }
  }

  removeInd <- unique(sort(removeDiffInd))
  keepInd <- setdiff(1:nrow(sv), removeInd)
  return(sv[keepInd])
}

selectDupSVtoRemove <- function(sv, dupInd1){
  removeInd <- NULL
  if (sv[dupInd1, sum(FILTER == "PASS")] >= 1 &&
      sv[dupInd1, sum(FILTER == "PASS")] < length(dupInd1)){
      removeInd <- dupInd1[sv[dupInd1, FILTER != "PASS"]] # remove not pass
  }else if (sv[dupInd1, sum(FILTER == "PASS")] == 0){ #none are PASS
    if (sv[dupInd1, length(unique(BXOL))] > 1){ # more than one BXOL value
      minADind <- sv[dupInd1, which(BXOL+BXSupport < max(BXOL+BXSupport))]
      rmInd <- dupInd1[minADind] # only keep highest support AD
    }else{ # if same AD, then keep max BXOL
      #rmInd <- dupInd1[sv[dupInd1, which(BXOL < max(BXOL))]]
      rmInd <- dupInd1 # remove all
    }
    removeInd <- rmInd
  }
  return(removeInd)
}

removeDupInterChrSV <- function(sv){
  ind <- sv[SPAN == -1, which = TRUE]
  dupInd <- which(duplicated(sv[ind, .(chromosome_1, start_1, chromosome_2, start_2)]))
  allInd <- 1:nrow(sv)
  sv <- sv[setdiff(allInd, ind[dupInd])]
  return(sv)
}

## removes 2nd instance of SVs with identical coordinates and orientation
removeIdenticalSV <- function(sv){
  dupInd <- which(duplicated(sv[, .(chromosome_1, start_1, chromosome_2, start_2, orient_1, orient_2)]))
  allInd <- 1:nrow(sv)
  sv <- sv[setdiff(allInd, dupInd)]
  return(copy(sv))
}

sortBkptPairOrder <- function(sv.unsort, chrs = c(1:22, "X", "Y")){
  genomeStyle <- seqlevelsStyle(sv.unsort$chromosome_1)[1]
  chrs <- as.character(chrs)
  seqlevelsStyle(chrs) <- genomeStyle
  ## interchr
  sv <- copy(sv.unsort)
  sv[, SV.id := 1:nrow(sv)]
  sv[, chromosome_1 := as.character(chromosome_1)]
  sv[, chromosome_2 := as.character(chromosome_2)]
  interChrInd <- sv[, which(chrs==chromosome_1) > which(chrs==chromosome_2), by=SV.id]$V1
  intraChrInd <- sv[, chromosome_1 == chromosome_2 & start_1 > start_2]
  ind <- which(interChrInd | intraChrInd)
  message("Re-assigning order for ", length(ind), " breakpoints.")
  #b1.colnames <- c("chromosome_1", "start_1", "alt_1", "orient_1", grep("\\.1", colnames(sv), value=T))
  #b2.colnames <- c("chromosome_2", "start_2", "alt_2", "orient_2", grep("\\.2", colnames(sv), value=T))
  b1.colnames <- grep("\\.1|\\_1", colnames(sv), value=T)
  b2.colnames <- grep("\\.2|\\_2", colnames(sv), value=T)
  sv.tmp <- cbind(sv[ind, b1.colnames, with=F])
  sv[ind, (b1.colnames) := sv[ind, b2.colnames, with=F]]
  sv[ind, (b2.colnames) := sv.tmp]
  return(copy(sv))
}

##########################################
##### remove 2nd mate ####### INCOMPLETE
##########################################
# pre-condition: MATEID info field with format nnnn:1 and nnnn:2
remove2ndMateCollapsedVCF <- function(vcf){
	indMATE <- grepl(":2", info(vcf)$MATEID)
	mate2 <- rowRanges(vcf)$ALT[indMATE]
	
	return(vcf)
}

##########################################
### Add new info field to VCF object #####
##########################################
addInfoFieldtoCollapsedVCF <- function(vcf, field.name, field.desc = "", 
															field.type = "Integer", field.defaultValue = NA, Number = 1){
	headerNames <- c(rownames(info(header(vcf))), field.name)
	suppressWarnings(info(header(vcf)) <- rbind(info(header(vcf)), 
									DataFrame(Number=Number, Type=field.type, Description=field.desc)))
	rownames(info(header(vcf))) <- headerNames
	newCol <- DataFrame(NA)
	colnames(newCol) <- field.name
	info(vcf) <- cbind(info(vcf), newCol)
	return(vcf)
}

########################
## Rsamtools BAM manipulations #
#########################
convertBamListToGRanges <- function(x){
	df <- as(as(x, "DataFrame"), "data.frame")
	seqnames <- as.character(df$rname)
	startPos <- df$pos
	endPos <- df$pos + df$qwidth
	df <- cbind(seqnames = df$rname, start = startPos, end = endPos, df)
	gr <- as(df, "GRanges")
	return(gr)	
}

### pre-condition: requires qname field
excludePartialMateMapping <- function(x, region, minMapQ = 60){
	## find reads that are fully contained "within" region of interest ##
	hits <- findOverlaps(query = x, subject = region, type = "within")
	x <- x[queryHits(hits), ]
	## remove reads less than minMapQ
	ind <- x$mapq >= minMapQ
	x <- x[ind, ]
	## remove reads if only present in region as a single mate (ie. other mate is aligned outside)
	ind <- duplicated(x$qname, fromLast = TRUE) | duplicated(x$qname, fromLast = FALSE)
	x <- x[ind, ]
	return(x)
}

## convert sv data.table to GenomicRanges
convert2GRanges <- function(sv, buffer = 0){
  sv.gr.1 <- copy(sv)
  setnames(sv.gr.1, c("chromosome_1", "start_1"), c("chr", "start"))
  sv.gr.1[, end := start + buffer]; sv.gr.1[, start := start - buffer]
  sv.gr.1 <- as(sv.gr.1, "GRanges")
  sv.gr.2 <- copy(sv)
  setnames(sv.gr.2, c("chromosome_2", "start_2"), c("chr", "start"))
  sv.gr.2[, end := start + buffer]; sv.gr.2[, start := start - buffer]
  sv.gr.2 <- as(sv.gr.2, "GRanges")
  return(list(sv.1 = sv.gr.1, sv.2 = sv.gr.2))
}

###########################################
### LOAD VCF FILES AS DATA.TABLE BY CHR ###
###########################################
loadVCFtoDataTableByChromosome <- function(svFiles, chrs = as.character(c(1:22, "X")), genomeStyle = "NCBI", genomeBuild = "hg19", minSPAN = 10000, minSPANBX = 10000, minBXOL = 0, maxBXOL = Inf, maxBXCount = 10000, minBSDsupport = 5, dupSV.bpDiff = 1000){
  seqlevelsStyle(chrs) <- genomeStyle
  if (length(svFiles) > 0){ ## sv vcf file exists ##
      vcf <- tryCatch({
        message("loading vcf ", svFiles)
        readVcf(svFiles, genome = genomeBuild)
      }, error = function(x){ 
        return(NA) 
      })
  }
  sv <- getSVfromCollapsedVCF(vcf, chrs = chrs, genomeStyle = genomeStyle)
  ## filter by BXOL support ##
  message("Filtering SVABA calls...")
  sv <- filterSVABA(sv, minSPAN=minSPAN, dupSV.bpDiff=dupSV.bpDiff)
   return(sv=sv)
}

filterSVABA <- function(x, minSPAN = 10000, dupSV.bpDiff = 1000){
  sv <- copy(x)
  ## filter by BXOL support ##
  sv$support <- NA
  indSnowman <- (sv$FILTER == "PASS")
  if (!is.null(sv$BXOL)){ ## if the barcode support is present
    indBarcode <- (sv$SPAN >= minSPANBX | sv$SPAN == -1) &
        ((sv$BXC.1 <= maxBXCount & sv$BXC.2 <= maxBXCount) | #both less than max count
        (sv$BXC.1 >= maxBXCount & sv$BXC.2 >= maxBXCount)) & #both higher than max count
        (sv$BXOL >= minBXOL & sv$BXOL <= maxBXOL) 
    sv$support[indBarcode] <- "BX"
  }else{
    indBarcode <- !logical(nrow(sv))
  }
  if (length(which(indSnowman | indBarcode)) == 0){
    stop("No SV events SVABA pass or bxol >= ", minBXOL)
  }
  sv$support[indSnowman] <- "SVABA"
  sv <- sv[which(indSnowman | indBarcode), ] 
  # concatenate all chromosomes
  #vcfBX <- vcf[which(indSnowman | indBarcode)]
  #if (j == 1){ vcfSampleBX <- vcfBX }
  #else{ vcfSampleBX <- rbind(vcfSampleBX, vcfBX); values(vcfSampleBX) <- values(vcfSampleBX)[1:5] }
  
  sv <- removeDupInterChrSV(sv)
  #sv <- removeDupNonPassSV(sv, bp.diff = dupSV.bpDiff)
  return(copy(sv))
}

setGenomeStyle <- function(x, genomeStyle = "NCBI", species = "Homo_sapiens"){
	#chrs <- genomeStyles(species)[c("NCBI","UCSC")]
	if (!genomeStyle %in% seqlevelsStyle(as.character(x))){
    	x <- suppressWarnings(mapSeqlevels(as.character(x), 
    					genomeStyle, drop = FALSE)[1,])
    }
    
    autoSexMChr <- extractSeqlevelsByGroup(species = species, 
    				style = genomeStyle, group = "all")
    x <- x[x %in% autoSexMChr]
    return(x)
}


###########################################
### LOAD BPS FILES AS DATA.TABLE BY CHR ###
###########################################
loadBPStoDataTableByChromosome <- function(bpsFile, tumor.id, chrs = c(1:22, "X"), 
    minLength = 1, dupSV.bpDiff = 1000){
  sv <- fread(paste0("gunzip -c ", bpsFile))
  sv <- sv[chr1 %in% chrs & chr2 %in% chrs]
  ## filter vcf by span length ##
  #indSPAN <- sv[, span >= minLength | span == -1 | confidence == "PASS"]
  ## filter vcf by FILTER == DUPREADS and LOCALMATCH
  #indFILTER <- sv[, !confidence %in% c("DUPREADS", "LOCALMATCH")]
  ## keep only 1st mate of each pair by looking at mate id of 2nd mate ##
  #ind <- indSPAN & indFILTER# & indMATE
  #message("Filtered ", sum(!ind), " rows - ", sum(ind), " remaining")
  #sv <- sv[ind]
  ## rename columns to match VCF
  setnames(sv, c("chr1","chr2","pos1","pos2","span","confidence"),
    c("chromosome_1","chromosome_2","start_1","start_2","SPAN","FILTER"))
  ## orientation: fwd = -, rev = +
  sv[, orient_1 := "fwd"]; sv[strand1 == "+", orient_1 := "rev"]
  sv[, orient_2 := "fwd"]; sv[strand2 == "+", orient_2 := "rev"]
  ## extract support counts from genotype column of tumor
  tumCol <- grepl(tumor.id,colnames(sv))
  geno <- do.call(rbind, strsplit(as.data.frame(sv[, tumCol, with=F])[,1], ":"))
  sv[, DR := as.integer(geno[, 7])]
  sv[, SR := as.integer(geno[, 6])]
  sv[, AD := as.integer(geno[, 2])]
  sv[, DP := as.integer(geno[, 3])]
  # remove duplicate events due to parallelizing analysis by chr
  #sv <- removeDupInterChrSV(sv)
  # remove unfiltered events that are close in coordinates to PASS events
  #sv <- removeDupNonPassSV(sv, bp.diff = dupSV.bpDiff)
  return(copy(sv))
}

annotVCFbyOverlap <- function(sv1, sv2, annotCol, buffer = 0){
  sv1.gr.1 <- copy(sv1)
  setnames(sv1.gr.1, c("chromosome_1","start_1"), c("chr","start"))
  sv1.gr.1[, end := start + buffer]; sv1.gr.1[, start := start - buffer]
  sv1.gr.1 <- as(sv1.gr.1, "GRanges")
  sv1.gr.2 <- copy(sv1)
  setnames(sv1.gr.2, c("chromosome_2","start_2"), c("chr","start"))
  sv1.gr.2[, end := start + buffer]; sv1.gr.2[, start := start - buffer]
  sv1.gr.2 <- as(sv1.gr.2, "GRanges")
  sv2.gr.1 <- copy(sv2)
  setnames(sv2.gr.1, c("chromosome_1","start_1"), c("chr","start"))
  sv2.gr.1[, end := start + buffer]; sv2.gr.1[, start := start - buffer]
  sv2.gr.1 <- as(sv2.gr.1, "GRanges")
  sv2.gr.2 <- copy(sv2)
  setnames(sv2.gr.2, c("chromosome_2","start_2"), c("chr","start"))
  sv2.gr.2[, end := start + buffer]; sv2.gr.2[, start := start - buffer]
  sv2.gr.2 <- as(sv2.gr.2, "GRanges")

  hits1 <- findOverlaps(query=sv1.gr.1, subject=sv2.gr.1)
  hits2 <- findOverlaps(query=sv1.gr.2, subject=sv2.gr.2)
  ind <- rbind(data.frame(hits1), data.frame(hits2))
  ind <- ind[duplicated(ind), ]
  #sv[ind$queryHits, eval(annotCol) := sv2[ind$subjectHits, get(annotCol)]]
  return(list(ind = ind$queryHits, values = sv2[ind$subjectHits, get(annotCol)]))
}

#########################################
## fit BXOL by length from SNOWMAN SV ###
#########################################
computeBXOLbinomialTest <- function(svSample, minBXOL=3, minSPAN=10000, minSPANBX=10000,
      se.level=0.95, loess.span=0.1, filter.quantile=0.95){
  bxol.fit.intraChr <- fitBXOLbyLength(svSample, interchr = FALSE, minBXOL = minBXOL, 
      minLength = minSPAN, minBXLength = minSPANBX, se.level = se.level, loess.span = loess.span, 
      filter.quantile = filter.quantile)
  bxol.fit.interChr <- fitBXOLbyLength(svSample, interchr = TRUE, minBXOL = minBXOL, 
      se.level = se.level, loess.span = loess.span, filter.quantile = filter.quantile)
  sv <- bxol.fit.intraChr$sv
  # interchr - use median BX.frac.min of SVABA interchr events #
  sv[SPAN == -1, fit.estimate := bxol.fit.interChr$sv[SPAN==-1 & support=="SVABA" & BXOL >= minBXOL, max(median(BX.frac.min), 0)]]
  sv[SPAN == -1, fit.estimate.se.low := fit.estimate]
  
  ## compute p-value for each breakpoint in the pair ##
  # use binomial exact test #
  sv[support=="BX" & !is.na(BXOL) & BXOL >= minBXOL & !is.na(fit.estimate), 
    BXOL.pval.1 := pbinom(q=BXOL-1, size=BXC.1, prob=fit.estimate, lower.tail=FALSE)]
  sv[support=="BX" & !is.na(BXOL) & BXOL >= minBXOL & !is.na(fit.estimate), 
    BXOL.pval.2 := pbinom(q=BXOL-1, size=BXC.2, prob=fit.estimate, lower.tail=FALSE)]
  return(list(sv=copy(sv), gp=bxol.fit.intraChr$gp))
}

###########################################
########## FIT BXOL BY LENGTH #############
###########################################
fitBXOLbyLength <- function(svSample, interchr = FALSE, minBXOL = 3, minLength = 0, minBXLength = 0,
    se.level = 0.95, loess.span = 0.3, minNtoFit = 10, filter.quantile = 0.95){
  if (length(svSample$BXC.1) == 0 || length(svSample$BXC.2) == 0 || length(svSample$BXOL) == 0){
    stop("Missing 1 of columns: BXOL or BXC.1 or BXC.2")
  }
  sv <- copy(svSample)
  sv[, BX.frac.1 := BXOL / BXC.1]
  sv[, BX.frac.2 := BXOL / BXC.2]
  sv[, BX.frac.min := pmin(BX.frac.1, BX.frac.2)]
  sv[, SPAN.log10 := log10(SPAN)]
  if (!interchr){ # intra chr #
    maxLength <- max(sv[SPAN > minSPAN & support == "SVABA", SPAN])
  }else{  # inter chr #
    minLength <- -2
    minBXLength <- -2
    maxLength <- 0
  }
  maxBXOL <- quantile(sv[support == "SVABA" & SPAN > minLength & SPAN < maxLength & BXOL >= minBXOL, BXOL], filter.quantile)
  # fit to only specified data #
  ind <- sv[support == "SVABA" & SPAN > minLength & SPAN < maxLength & 
          BXOL >= minBXOL & BXOL <= maxBXOL, which=TRUE]
  fit.loess <- NULL; gp <- NULL
  if (!interchr){
    
    fit.loess <- tryCatch({
      sv[ind, loess(BX.frac.min ~ SPAN.log10, span = loess.span)]
    }, error = function(x){ 
      message("BXOL loess fit failed - using median PASS BX.FRAC.MIN instead.")
      return(NULL) 
    })
    if (!is.null(fit.loess) && length(ind) >= minNtoFit && sum(is.na(fit.loess$fitted)) < fit.loess$n){
      fit.pred <- predict(fit.loess, se=TRUE)
      sv[ind, fit := fit.pred$fit]
      sv[ind, fit.df := fit.pred$df]
      sv[ind, fit.se.high := fit + qt(se.level, fit.df) * fit.pred$se.fit]
      sv[ind, fit.se.low := fit - qt(se.level, fit.df) * fit.pred$se.fit]
     
      # predict barcode overlap count for non-SVABA SVs based on fit from SVABA SVs 
      fit.est <- predict(fit.loess, sv$SPAN.log10, se=TRUE)
      sv[, fit.estimate := fit.est$fit]
      sv[, fit.estimate.df := fit.est$df]
      sv[, fit.estimate.se.high := fit.estimate + qt(se.level, fit.estimate.df) * fit.est$se.fit]
      sv[, fit.estimate.se.low := pmax(fit.estimate - qt(se.level, fit.estimate.df) * fit.est$se.fit, 0)]
      maxLength <- sv[ind, max(SPAN, na.rm=T)] + 1 # reassign
      # for events larger than maxLength of SVABA, assign original BX.frac.min
      sv[SPAN >= maxLength & support == "BX", fit.estimate := as.numeric(BX.frac.min)]
      sv[SPAN >= maxLength & support == "BX", fit.estimate.se.high := as.numeric(BX.frac.min)]
      sv[SPAN >= maxLength & support == "BX", fit.estimate.se.low := pmax(as.numeric(BX.frac.min), 0)]
    }
    if (is.null(fit.loess) || length(ind) < minNtoFit || sum(is.na(fit.loess$fitted)) >= fit.loess$n ||
        sum(is.na(fit.est$fit)) == nrow(sv)){ # fit did not work - not enough SVABA PASS data points
      # if no fit value assigned, then use sample median BX.FRAC.MIN for each SV
      bxInd <- sv[SPAN > minBXLength & SPAN < maxLength & BXOL >= minBXOL & support == "BX", which=T]
      medianSVABA.bx.frac.min <- sv[bxInd, max(median(BX.frac.min), 0)]
      sv[bxInd, fit.estimate := medianSVABA.bx.frac.min]
      sv[bxInd, fit.estimate.se.high := medianSVABA.bx.frac.min]
      sv[bxInd, fit.estimate.se.low := medianSVABA.bx.frac.min]
    }
    # plot data use for fitting and loess fit #
    gp <- ggplot(sv[SPAN < maxLength & BXOL >= minBXOL & BXOL <= maxBXOL &
          ((SPAN > minLength & support=="SVABA") | (SPAN > minBXLength & support=="BX"))][order(support)], 
          aes(x=SPAN, y=BX.frac.min, color=support, shape=support, alpha=support)) + 
        geom_point(size=1) + scale_alpha_manual(values=c(0.4, 0.8)) +
        scale_x_continuous(trans="log10") + 
        scale_color_manual(values=c("grey","black")) + scale_shape_manual(values=c(19,19)) +
        ylab("Barcode Overlap Fraction") + xlab("Breakpoint Pair Distance (bp)") + theme_bw() + 
        coord_cartesian(ylim=c(-0.1, 1.1)) + 
        theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
            axis.title=element_text(size=20), 
            plot.title=element_text(size=24, face="bold")) 
        #geom_smooth(method="loess", se=TRUE, span=0.3, size=1.5) 
    if ("fit.se.low" %in% colnames(sv)){
      gp <- gp + 
        geom_ribbon(data=sv[ind], aes(ymin=fit.se.low, ymax=fit.se.high), alpha=0.4, color=NA, fill="red") +
        geom_line(data=sv[ind], aes(y=fit), color="blue", size=1) + guides(alpha=FALSE) 
        #geom_line(data=sv[ind], aes(y=fit.se.high), linetype="dashed", color="grey", size=0.75) +
        #geom_line(data=sv[ind], aes(y=fit.se.low), linetype="dashed", color="grey", size=0.75) +      
    }
  
  }
  return(list(fit.loess=fit.loess, gp=gp, sv=copy(sv), max.length=maxLength, max.BXOL=maxBXOL))
}

###############################################
############# FIND SV EVENT TYPE ##############
###############################################
getSVType <- function(sv, minColSPAN = 10000, minTrans = 1e7){#, maxInvSPAN = 1e6, maxFBISPAN = 30000){
	sv <- copy(sv)
	#sv[SPAN > minColSPAN, type := NA]
	sv[orient_1=="rev" & orient_2=="fwd" & SPAN > minColSPAN, type := "Deletion"]
	sv[orient_1=="fwd" & orient_2=="rev" & SPAN > minColSPAN, type := "Duplication"]
	sv[(orient_1==orient_2) & SPAN > 0, type := "Inversion"]
	#sv[orient_1=="fwd" & orient_2=="fwd" & SPAN < maxFBISPAN & SPAN > 0, type := "FoldBackInv-T"]
	#sv[orient_1=="rev" & orient_2=="rev" & SPAN < maxFBISPAN & SPAN > 0, type := "FoldBackInv-H"]
	sv[SPAN == -1, type := "InterChr"]
	sv[SPAN > minTrans, type := "Translocation"]
	return(sv$type)
}


#################################################
############## OVERLAP SEGS AND SV ##############
#################################################
getSegSVoverlap <- function(segs, sv, event.cn, buffer = 1e5, interChr = FALSE){
	segs <- segs[Corrected_Call %in% event.cn]
	if (interChr){
	  sv <- sv[SPAN == -1]
	}else{
  	sv <- sv[SPAN > 0]
  }
	seg.sv <- NULL
	if (nrow(segs) > 0 && nrow(sv) > 0){
		segs_1 <- copy(segs)
		segs_1[, c("Start", "End") := data.frame(Start - buffer, Start + buffer)]
		segs_1.gr <- as(segs_1, "GRanges")
		segs_2 <- copy(segs)
		segs_2[, c("Start", "End") := data.frame(End - buffer, End + buffer)]
		segs_2.gr <- as(segs_2, "GRanges")

		sv_1 <- copy(sv)
		sv_1[, end := start_1]
		setnames(sv_1, c("chromosome_1","start_1"), c("chromosome", "start")) 
		sv_1.gr <- as(sv_1, "GRanges")
		sv_2 <- copy(sv)
		sv_2[, end := start_2]
		setnames(sv_2, c("chromosome_2","start_2"), c("chromosome", "start")) 
		sv_2.gr <- as(sv_2, "GRanges")

    #hits.a.b - a is sv breakpoint 1/2 and b is seg start/end
		hits1.1 <- findOverlaps(subject = segs_1.gr, query = sv_1.gr)
		hits1.2 <- findOverlaps(subject = segs_2.gr, query = sv_1.gr)
		hits2.1 <- findOverlaps(subject = segs_1.gr, query = sv_2.gr)
		hits2.2 <- findOverlaps(subject = segs_2.gr, query = sv_2.gr)
		#segs$bkpt1.1 <- FALSE; segs$bkpt1.2 <- FALSE
		#segs$bkpt2.1 <- FALSE; segs$bkpt2.2 <- FALSE

		#segs[subjectHits(hits1.1), bkpt1.1 := TRUE]; segs[subjectHits(hits1.2), bkpt1.2 := TRUE]
		#segs[subjectHits(hits2.1), bkpt2.1 := TRUE]; segs[subjectHits(hits2.2), bkpt2.2 := TRUE]
		#overlapInd <- segs.gain[, .I[bkpt1==TRUE | bkpt2==TRUE]]
		overlapInd <- unique(cbind(
		  c(subjectHits(hits1.1), subjectHits(hits1.2), subjectHits(hits2.1), subjectHits(hits2.2)), 
		  c(queryHits(hits1.1), queryHits(hits1.2), queryHits(hits2.1), queryHits(hits2.2))))
		seg.sv <- cbind(segs[overlapInd[, 1]], sv[overlapInd[, 2]])
    seg.sv[chromosome_1==Chromosome, bkpt1.1 := abs(start_1 - Start)]
    seg.sv[chromosome_1==Chromosome, bkpt1.2 := abs(start_1 - End)]
    seg.sv[chromosome_2==Chromosome, bkpt2.1 := abs(start_2 - Start)]
    seg.sv[chromosome_2==Chromosome, bkpt2.2 := abs(start_2 - End)]
			
		if (!interChr){
      if (nrow(seg.sv) > 0){ 
        # find overlap length #
        segs.gr <- as(seg.sv[, .(Chromosome, Start, End)], "GRanges")
        sv.gr <- seg.sv[, .(chromosome_1, start_1, start_2)]
        setnames(sv.gr, c("chromosome", "start", "end"))
        sv.gr <- as(sv.gr, "GRanges")
        hits0 <- findOverlaps(subject = segs.gr, query = sv.gr)
        overlap.len <- width(overlapsRanges(hits0, query = ranges(sv.gr), subject = ranges(segs.gr)))
        names(overlap.len) <- queryHits(hits0)
        sameRowInd <- which(queryHits(hits0) == subjectHits(hits0))
        seg.sv[queryHits(hits0)[sameRowInd], overlap.length := overlap.len[sameRowInd]]
        setkey(seg.sv, SEG.id)
        #seg.sv[, keep.sv := {
        #	ind = TRUE
        #	ind2 = rep(TRUE, .N)
        #	if (.N > 1){  # when more than breakpoint overlaps, pick the one with PASS or higher evidence
        #		ind = FILTER == "PASS" 
            #if (sum(ind) == 0){
            #	ind = (BXOL + AD) == pmax(BXOL + AD)
            #}
        #	}
        #	ind & ind2
        #}, by=c("Start","End")]
        seg.sv[, sv.overlap.frac := overlap.length / (SPAN + 1)]
      }else{
        seg.sv <- NULL
      }
    }
	}
	return(copy(seg.sv))
}

#########################################################
########## GET TANDEMDUP - FLANKING CN IS LOWER #########
#########################################################
# true if left and right flanking CN is lower 
# segs must have Chromosome and Sample columns
getTDcnFlank <- function(segs, maxDupLength = 10e7){
  cn <- copy(segs)
	cn[, flankInd := {
	    length = End - Start + 1           
      lengthInd = length < maxDupLength
      leftDiff = Copy_Number - c(NA, Copy_Number[-.N])
      leftInd = leftDiff > 0 #& leftLen
      rightDiff = Copy_Number - c(Copy_Number[-1], NA)
      rightInd = rightDiff > 0 #& rightLen
      flankDiff = abs(c(NA, Copy_Number[-.N]) - c(Copy_Number[-1], NA))
      flankInd = flankDiff <= 1#MAD/2
      leftInd & rightInd & flankInd & lengthInd
      }, by=c("Chromosome", "Sample")]
  return(cn[, flankInd])
}

############################################
########## ANNOTATION SV WITH CN ###########
############################################
annotateSVwithCN <- function(sv, cn, cnColToAnnotate = "Copy_Number", offset = 0){
  cn.gr <- as(cn, "GRanges")
  sv.1 <- copy(sv)
  setnames(sv.1, c("chromosome_1"), c("Chr"))
  sv.1[, start := start_1 + offset]
  sv.1[, end := start]
  sv.1 <- as(sv.1, "GRanges")
  sv.2 <- copy(sv)
  setnames(sv.2, c("chromosome_2"), c("Chr"))
  sv.2[, start := start_2 + offset]
  sv.2[, end := start]
  sv.2 <- as(sv.2, "GRanges")

  hits1 <- findOverlaps(query = sv.1, subject = cn.gr)
  bkpt1 <- cn[subjectHits(hits1), get(cnColToAnnotate)]
  hits2 <- findOverlaps(query = sv.2, subject = cn.gr)
  bkpt2 <- cn[subjectHits(hits2), get(cnColToAnnotate)]
  
  return(list(annot1 = bkpt1, ind1 = queryHits(hits1), annot2 = bkpt2, ind2 = queryHits(hits2)))
}

#######################################################
########## ANNOTATION SV WITH CN FROM NEAREST BIN ###########
#######################################################
## modified to work for combine_TITAN_ICHOR/titan_ichor_cn.txt
## returns datatable for prev and next (flanking) copy number
## returns datatable with indices used in "cn" input
annotateSVwithFlankCN <- function(sv, cn, cnColToAnnotate = "Copy_Number", buffer = 1e6){
	x <- copy(sv); y <- copy(cn)
	z <- data.table()
	indCN <- data.table()
	#y <- na.omit(y)
	for (i in 1:nrow(sv)){ # for each sv
    bin1Ind <- x[i, chromosome_1] == y[, Chr] & 
                x[i, start_1] >= (y[, Start]) & 
                x[i, start_1] <= (y[, End] + buffer)
    bin2Ind <- x[i, chromosome_2] == y[, Chr] & 
                x[i, start_2] >= (y[, Start]) & 
                x[i, start_2] <= (y[, End] + buffer)
    # look at closest left and right points, find minimum distance to points, assign copy of min point
    ind1.prev <- NA; ind1.next <- NA; ind2.prev <- NA; ind2.next <- NA;
    y1.prev <- NA; y1.next <- NA; y2.prev <- NA; y2.next <- NA;
    dist1.prev <- NA; dist1.next <- NA; dist2.prev <- NA; dist2.next <- NA;
    if (sum(bin1Ind) > 0){
      ind <- tail(which(bin1Ind), 1)
      ind1.prev <- max(ind - 1, 1)
      ind1.next <- min(ind + 1, length(bin1Ind)) 
      y1.prev <- y[ind1.prev, get(cnColToAnnotate)]
      dist1.prev <- abs(x[i, start_1] - y[ind1.prev, End])
      y1.next <- y[ind1.next, get(cnColToAnnotate)]
      dist1.next <- abs(x[i, start_1] - cn[ind1.next, Start])
      if (y[ind, Chr] != y[ind1.prev, Chr]){
        ind1.prev <- NA; y1.prev <- NA; dist1.prev <- NA;
      }
      if (y[ind, Chr] != y[ind1.next, Chr]){
        ind1.next <- NA; y1.next <- NA; dist1.next <- NA;
      }
    }
    if (sum(bin2Ind) > 0){
      ind <- tail(which(bin2Ind), 1)
      ind2.prev <- max(ind - 1, 1)
      ind2.next <- min(ind + 1, length(bin1Ind))
      y2.prev <- y[ind2.prev, get(cnColToAnnotate)]
      dist2.prev <- abs(x[i, start_2] - y[ind2.prev, End])
      y2.next <- y[ind2.next, get(cnColToAnnotate)]
      dist2.next <- abs(x[i, start_2] - y[ind2.next, Start])
      if (y[ind, Chr] != y[ind2.prev, Chr]){
        ind2.prev <- NA; y2.prev <- NA; dist2.prev <- NA;
      }
      if (y[ind, Chr] != y[ind2.next, Chr]){
        ind2.next <- NA; y2.next <- NA; dist2.next <- NA;
      }
	  } 
	  indCN <- rbind(indCN, data.table(ind1.prev, ind1.next, ind2.prev, ind2.next))
	  z <- rbind(z, data.table(y1.prev, dist1.prev, y1.next, dist1.next, y2.prev, dist2.prev, y2.next, dist2.next))
	}
	colnames(z) <- paste0(cnColToAnnotate, c("_1_prev", "_1_prev_dist", "_1_next", "_1_next_dist", 
	  "_2_prev", "_2_prev_dist",  "_2_next", "_2_next_dist"))
	colnames(indCN) <- paste0(cnColToAnnotate, c("_1_prev", "_1_next", "_2_prev", "_2_next"))
	return(list(annot = z, ind= indCN))
}

#######################################################
########## ANNOTATION BETWEEN SV BREAKPOINTS ###########
#######################################################
# for each pair of SV breakpoints, what is the annotation between them?
annotateSVbetweenBkptsWithCN <- function(sv, cn, segs, buffer = 1e5, 
    cnColToAnnotate = "Copy_Number", fun = "mean"){
  ind <- sv[chromosome_1 == chromosome_2 & SPAN > 0, which = T]
  sv.intraChr <- copy(sv[ind])
  sv.gr <- copy(sv.intraChr)
  setnames(sv.gr, c("chromosome_1", "start_1", "start_2"), c("Chr", "Start", "End"))
  sv.gr <- as(sv.gr, "GRanges")
  cn.gr <- as(cn, "GRanges")
  segs.gr <- as(segs, "GRanges")
  hits.cn <- findOverlaps(subject = sv.gr, query = cn.gr, type = "within")
  annot.cn <- ldply(as.list(t(hits.cn)), .fun=function(x){ 
    cn <- cn[x, do.call(fun, list(get(cnColToAnnotate), na.rm=T))] 
    data.table(N=length(x), cn)
  })
  annot.cn <- data.table(annot.cn)
  hits.seg <- findOverlaps(query = sv.gr, subject = segs.gr, type = "any")  
  hits.seg.l <- as.list(hits.seg)
  annot.seg <- ldply(1:length(hits.seg.l), .fun=function(i){ 
    x <- hits.seg.l[[i]]
    # get overlap lengths - only look at segments > buffer
     ovLen <- width(overlapsRanges(hits.seg[queryHits(hits.seg)==i], query=ranges(sv.gr), subject=ranges(segs.gr)))
    ov.ind <- ovLen > buffer
    x <- x[ov.ind]; ovLen <- ovLen[ov.ind]    
    id <- segs[x, paste0(SEG.id, collapse=",")] 
    cn <- segs[x, paste0(get(cnColToAnnotate), collapse=",")]
    data.table(NumSeg = length(x), cn, SEG.id=id)
  })
  annot.seg <- data.table(annot.seg)
  return(list(annot.cn=annot.cn, annot.seg=annot.seg, ind=ind))
}


###########################################################
########## FIND CN BETWEEN ADJACENT BREAKPOINTS ###########
###########################################################
### INCOMPLETE ###
findDeletionBridges <- function(sv, cn, overlapFrac = 0.5){
  cn.gr <- as(cn, "GRanges")
  sv.1.next <- copy(sv)
  setnames(sv.1.next, c("chromosome_1", "start_1"), c("Chr","start"))
  sv.1.next[, end := {
    ind = Chr == c(Chr[-1], NA)
    adj.chr = c(start)
  }] ## set next SV start_1 as end
  
  
  sv.prev <- copy(sv)
  sv.prev[, end := c(NA, start_2[-.N])] ## set next SV start_1 as end
  setnames(sv.1, c("chromosome_1"), c("Chr"))
  sv.1[, start := start_1 + offset]
  sv.1[, end := start]
  sv.1 <- as(sv.1, "GRanges")
  sv.2 <- copy(sv)
  setnames(sv.2, c("chromosome_2"), c("Chr"))
  sv.2[, start := start_2 + offset]
  sv.2[, end := start]
  sv.2 <- as(sv.2, "GRanges")

}

###########################################
############ PLOTTING FUNCTION ###########
###########################################
### TO-DO: 
# 1) Filter by support
# 2) Strand direction of breakpoint
plotRearrangementArcs <- function (vcf, hmmcopy, ploidy = NULL, interchr=TRUE, xlim=NULL, segIn=NULL, support=1, chr=NULL, chrLens=NULL, minSPAN = 10, buffer=1e6, centreLine=NULL, arcHeight=4, lty=1, lcol = "black", svTypeCol=FALSE, arr.col = "black", arr.pos = 1, endhead = FALSE, lwd = 1, orient="topbottom", include.inter.chr = FALSE){
	require(diagram)
	if (class(vcf)[1] == "CollapsedVCF"){
		sv <- getSVfromCollapsedVCF(vcf)
		sv$support <- NA
		indSnowman <- (sv$FILTER == "PASS")
		sv$support[indSnowman] <- "SVABA"
		if (!interchr){
		  sv <- sv[sv$SPAN > 0, ]
		}
		if (!is.null(sv$BXOL)){ ## if the barcode support is present
			indBarcode <- (sv$SPAN >= minSPAN) & # | sv$SPAN == -1) &
						(sv$BXOL >= minBXOL & sv$BXOL <= maxBXOL & (sv$BXOL+sv$AD) >= minBSDsupport)
			sv$support[indBarcode] <- "BX"
			sv <- sv[indSnowman | indBarcode, ]
		}else{
			indSnowman <- sv$SPAN >= minSPAN | sv$SPAN == -1
			sv <- sv[indSnowman, ]
		}		
	}else if (is.data.frame(sv)){
		sv <- vcf
	}
	sv <- sv[sv$SPAN >= minSPAN | sv$SPAN == -1, ]
	sv <- as.data.frame(sv)
	 #remove duplicate breakpoints
	#sv <- removeDupSV(sv)
	
		# adjust for ploidy #
	if (!is.null(ploidy)){
    hmmcopy[, "copy"] <- as.numeric(hmmcopy[, "copy"]) + log2(ploidy / 2)
  }
	
  if (!is.null(chr) && nrow(sv) > 0){
    dataByChr <- sv[sv[,"chromosome_1"]==as.character(chr) | 
                  sv[,"chromosome_2"]==as.character(chr),]
    if (nrow(dataByChr) > 0){
      #use only interchromosomal rearrangements
      dataByChr$intraChr <- dataByChr[,"chromosome_1"] == dataByChr[,"chromosome_2"]
      dataByChr$inWindow <- "NONE"
      if (!is.null(xlim)){ 
        ## only include SV if at least one of the breakpoints is inside window
        leftIn <- (dataByChr[, "chromosome_1"] == chr & 
                      dataByChr[, "start_1"] >= xlim[1] - buffer & 
                      dataByChr[, "start_1"] <= xlim[2] + buffer)
        rightIn <- (dataByChr[, "chromosome_2"] == chr & 
                      dataByChr[, "start_2"] >= xlim[1] - buffer & 
                      dataByChr[, "start_2"] <= xlim[2] + buffer)
        dataByChr$inWindow[leftIn | rightIn] <- "ONE"
        ## SV for both breakpoints contained within window
        dataByChr$inWindow[leftIn & rightIn] <- "BOTH"						
      }		
		#intraChrSV <- dataByChr[intrachrInd, ]
		#interChrSV <- dataByChr[interchrInd, ]
		
		#hmmcopy$midCoord <- (hmmcopy$end - hmmcopy$start + 1) / 2 + hmmcopy$start
		## plot intra-chromosomal SV ##
		
			for (i in 1:nrow(dataByChr)){
				if (dataByChr[i, "inWindow"] == "NONE"){
					next
				}
				#print(dataByChr[i, ])
				## get hmmcopy logR of nearest bin
				yLogR <- findNearestLogR(dataByChr[i, ], hmmcopy, buffer = buffer * 100)	
				yToUse <- yLogR			
				orient <- "top"
				start1 <- dataByChr[i,"start_1"]
        end1 <- dataByChr[i,"start_2"]
				#if (dataByChr[i, "intraChr"] && dataByChr[i, "inWindow"] == "BOTH"){ # intra chr event #
					## breakpoint type ##
          # TANDEM DUPLICATION #
#					if (dataByChr[i, "orient_1"] == "fwd" && dataByChr[i, "orient_2"] == "rev"){
#						start1 <- dataByChr[i,"start_2"]
#						end1 <- dataByChr[i,"start_1"]
#						yToUse <- yLogR[c(2,1)]
#						arr.col.sv <- "red"
#						sv.type <- "dup"
					## DELETION ##
#					}else if (dataByChr[i, "orient_1"] == "rev" && dataByChr[i, "orient_2"] == "fwd"){
#						start1 <- dataByChr[i,"start_1"]
#						end1 <- dataByChr[i,"start_2"]
#						yToUse <- yLogR
#						arr.col.sv <- "green"
#						sv.type <- "del"
					## INVERSION ## - will include 2 records, but will overlap in plot 
#					}else if ((dataByChr[i, "orient_1"] == "fwd" && dataByChr[i, "orient_2"] == "fwd") ||
#					          (dataByChr[i, "orient_1"] == "rev" && dataByChr[i, "orient_2"] == "rev")){
#					  start1 <- dataByChr[i,"start_1"]
#					  end1 <- dataByChr[i,"start_2"]
#					  yToUse <- yLogR[c(1,2)]
#					  start2 <- dataByChr[i,"start_2"]
#					  end2 <- dataByChr[i,"start_1"]
#					  yToUse2 <- yLogR[c(2,1)]
#					  arr.col.sv <- "blue"
#					  sv.type <- "inv"
#					}
#					if (dataByChr[i, "SPAN"] == -1){
#					  start1 <- dataByChr[i,"start_1"]
#						end1 <- dataByChr[i,"start_2"]
#						yToUse <- yLogR[c(1,2)]
#					  arr.col.sv <- "black"
#					  sv.type <- "interchr"
#					}
					if (!is.null(centreLine)){
					  yToUse <- c(centreLine, centreLine)
					  yToUse2 <- c(centreLine, centreLine)
					}
					curveHeight <- arcHeight - mean(yToUse)#/ abs(start1 - end1)
					if (is.null(arr.col)){
					  arr.col.toUse <- arr.col.sv
					}else {arr.col.toUse <- arr.col}
					if (is.null(lcol)){
				    lcol.toUse <- arr.col.sv
				  }else {lcol.toUse <- lcol}
					if (svTypeCol){
					  lcol.toUse <- dataByChr[i, "color"]
					}
					## if only 1 breakpoint in window ##
					if (dataByChr[i, "inWindow"] == "ONE"){
					## if inter-chromosomal ##
		        if (dataByChr[i, "chromosome_1"] != chr || dataByChr[i, "chromosome_2"] != chr){
		          orient <- "bottom"
		        }
		
		       ## for whole chromosome
					  if (dataByChr[i, "chromosome_1"] != chr){					    
              if (abs(dataByChr[i, "start_2"]-xlim[1]) < abs(dataByChr[i, "start_2"]-xlim[2])){
              # 2nd bkpt is closer to chr begin
                start1 <- par("usr")[1] - (par("usr")[2] * 1.5 - par("usr")[2])# left boundary  
              }else{
                start1 <- end1
                end1 <- par("usr")[2]*1.5 #right boundary  
                yToUse <- yLogR[c(2,1)]              
              }                   
            }else if (dataByChr[i, "chromosome_2"] != chr){
              if (abs(dataByChr[i, "start_1"]-xlim[1]) < abs(dataByChr[i, "start_1"]-xlim[2])){
              # 1st bkpt is closer to chr begin
                end1 <- start1
                start1 <- par("usr")[1] - (par("usr")[2] * 1.5 - par("usr")[2]) # left boundary 
                yToUse <- yLogR[c(2,1)]               
              }else{
                end1 <- par("usr")[2]*1.5 #right boundary      
              }  
            }
            ## for zoom within same chromosome
            if (dataByChr[i, "intraChr"]){
              if (start1 < xlim[1]){ 
                start1 <- par("usr")[1] - (par("usr")[2] * 1.5 - par("usr")[2]) # left boundary  
              }
              if (start1 > xlim[2]){
                start1 <- par("usr")[2] - (par("usr")[2] * 1.5 - par("usr")[2]) #right boundary
              }
              if (end1 < xlim[1]){
                end1 <- par("usr")[1]*1.5 # left boundary
              }
              if (end1 > xlim[2]){
                end1 <- par("usr")[2]*1.5 # left boundary               
              }        
              #if (sv.type == "inv") { end2 <- start1; start2 <- end1 }    
            }
            curveHeight <- par("usr")[3] - mean(yToUse)#* 2 / (end1 - start1)
					}
					
					curvedarrow.new(from = c(start1, yToUse[1]), endhead = endhead, orient = orient,
										to = c(end1, yToUse[2]), lty = lty, lwd = lwd, lcol=lcol.toUse, arr.col=arr.col.toUse, 
										arr.width = 0.1, arr.length=0.2,
										curve = curveHeight, arr.type = "triangle", arr.pos = arr.pos)
					#if (sv.type == "inv" && dataByChr[i, "inWindow"] == "BOTH"){ ## plot second arc for inversions
					#  curvedarrow.new(from = c(start2, yToUse2[1]), endhead = endhead, orient = orient,
					#					to = c(end2, yToUse2[2]), lty = lty, lwd = lwd, lcol=lcol.toUse, arr.col=arr.col.toUse, 
					#					arr.width = 0.1,arr.length=0.2,
					#					curve = curveHeight, arr.type = "triangle", arr.pos = arr.pos)
					#}
					
					#iArrows(start1, yToUse[1], end1, yToUse[2], h.lwd=2, sh.lwd=2, sh.col=lcol, curve=curveHeight,
					#				width=1, size=0)
					## plot 2nd breakpoint ##
					#if (dataByChr[i, "orient_2"] == "fwd"){
					#	start2 <- dataByChr[i,"start_1"]
					#	end2 <- dataByChr[i,"start_2"]
					#	yToUse <- yLogR
					#}else if (dataByChr[i, "orient_2"] == "rev"){
					#	start2 <- dataByChr[i,"start_2"]
					#	end2 <- dataByChr[i,"start_1"]
					#	yToUse <- yLogR[c(2,1)]
					#}
					#curvedarrow(from = c(start2, yToUse[1]), endhead = endhead,
					#					to = c(end2, yToUse[2]), lwd = 1, lcol=lcol, arr.col=arr.col,
					#					curve = curveHeight, arr.type = "triangle", arr.pos = 1 - arr.pos)
				#}else if (dataByChr[i, "inWindow"] == "ONE"){ # inter chr event #
				#	orient <- "top"
				#	if (dataByChr[i, "chromosome_1"] == dataByChr[i, "chromosome_2"]){
				#		start1 <- dataByChr[i,"start_1"]
				#		end1 <- dataByChr[i,"start_2"]
				#		orient <- "top"
				#	}else if (dataByChr[i, "chromosome_1"] != chr){
				#		start1 <- par("usr")[2]*1.05
				#		end1 <- dataByChr[i,"start_2"]
				#	}else{
				#		end1 <- par("usr")[2]*1.05
				#		start1 <- dataByChr[i,"start_2"]
				#	}
				#	if (!is.null(centreLine)){
				#	  yToUse <- c(centreLine,centreLine)
				#	}else{
				#	  yToUse <- yLogR
				#	}
				#	curveHeight <- arcHeight * 2 / (end1 - start1)
				#	curvedarrow.new(from = c(start1, yToUse[1]), endhead = endhead, orient = orient,
				#					to = c(end1, yToUse[2]), lwd = lwd, lcol=lcol, arr.col=arr.col, arr.lwd = 0.5,
				#					curve = -1 * curveHeight, arr.type = "triangle", arr.pos = 0.5)
				#}
			}
		}		
	}
	invisible()
}

## compute copy number using corrected log ratio ##
logRbasedCN <- function(x, purity, ploidyT, cn = 2){
	ct <- (2^x * (cn * (1 - purity) + purity * ploidyT * (cn / 2)) - cn * (1 - purity)) / purity
	ct <- sapply(ct, max, 1/2^6)
	return(ct)
}


## format of dataIn is output from /home/unix/gavinha/software/code/git/scripts/titan/analysis/combineTITAN-ichor.R
plotTitanIchorCNA <- function(dataIn, param = NULL, colName = "LogRatio", segs=NULL, chr=NULL, purity = NULL, ploidyT = NULL, geneAnnot=NULL, yrange=c(-4,6), yaxis = "logRatio", xlim=NULL, xaxt = "n", cex = 0.5, gene.cex = 0.5, plot.title = NULL, cnCol = NULL, spacing=4, cytoBand=T, alphaVal=1, main){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  alphaSubcloneVal <- ceiling(alphaVal / 2 * 255); class(alphaVal) = "hexmode"
  subcloneCol <- c("#00FF00")
  if (is.null(cnCol)){
		cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 26))
	}
	cnCol <- paste(cnCol,alphaVal,sep="")
  names(cnCol) <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
	cnCol <- c(cnCol, "HET"="#0000FF", "DLOH"="#006400", "NLOH"="#0000FF", "ALOH"="#FF0000", "ASCNA"="#FF0000", "BCNA"="#FF0000", "UBCNA"="#FF0000")
  # adjust for ploidy #
  normCN <- 2
  if (!is.null(ploidyT)){
    ploidyS <- purity * ploidyT + (1-purity) * normCN
   # dataIn[, colName] <- as.numeric(dataIn[, colName])# + log2(ploidyS / 2)
    
    if (!is.null(segs)){
      segs[, colName] <- segs[, colName] + log2(ploidyS / 2)
    }
  }
  
  
  if (!is.null(chr)){
    for (i in chr){
      dataByChr <- dataIn[dataIn[,"Chr"]==as.character(i),]
       ## set y axis labels as either integer or logR copy number
      #avgTumPloidy <- round(ploidyT)
      zero <- 0.25  
      cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      #ploidyToUse <- ploidyS
      if (i == "X"){
        normCN <- 1
        zero <- 0.5
        cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      }      
      if (yaxis == "integer"){
        y.ticks <- log2(cn)
        y.ticks[1] <- log2(zero)  
        yrange[1] <- y.ticks[1]    
        ylab <- "Copy Number"
        #dataByChr[, colName] <- log2(logRbasedCN(dataByChr[, colName], purity, ploidyT, cn=normCN))
        dataByChr[, colName] <- log2(dataByChr[, colName])
        if (!is.null(segs)){
      		segs[, colName] <- log2(segs[, colName])# + log2(ploidyS / 2)
    		}
        centreLine <- log2(normCN)
      }else{     
        cn <- (yrange[1]:yrange[2])
      	dataByChr[, colName] <- dataByChr[, colName] + log2(ploidyS / 2)
        y.ticks <- cn
        ylab <- "Copy Number (log2 ratio)"
        centreLine <- 0
      }

      #plot the data
      #if (outfile!=""){ pdf(outfile,width=10,height=6) }
      par(mar=c(spacing,8,4,2))
      #par(xpd=NA)
      coord <- (as.numeric(dataByChr[,"End"]) + as.numeric(dataByChr[,"Start"]))/2
      if (is.null(xlim)){
        xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"Start"]))
        xaxt <- "n"
      }
      if (is.null(plot.title)){
        plot.title <- paste("Chromosome ",i,sep="")
      }
      ## plot logR for bins ##
      plot(coord,as.numeric(dataByChr[, colName]),col=cnCol[dataByChr[,"event"]],
           pch=16, ylim=yrange, yaxt="n",
           xlim=xlim, xaxt = xaxt, xlab="",ylab=ylab,
           cex.lab=1.5,cex.axis=1.5, cex=cex,las=1)
      axis(2, at=y.ticks, labels=cn, las=2, cex.axis=1.5)
      title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)
      ## plot centre line ##
      lines(c(1,as.numeric(dataByChr[dim(dataByChr)[1],3])),rep(centreLine,2),type="l",col="grey",lwd=0.75)
      if (!is.null(segs)){
        segsByChr <- segs[segs[,"Chromosome"]==as.character(i),,drop=FALSE]
        #ind <- segsByChr$subclone.status == FALSE
        apply(segsByChr, 1, function(x){
          lines(x[c("Start","End")], rep(x[colName], 2), col = cnCol[x["event"]], lwd = 3)
          invisible()
        })
        #if (sum(!ind) > 0){
        #  apply(segsByChr[!ind, ], 1, function(x){
        #    lines(x[c("Start","End")], rep(x["Median_logR"], 2), col = subcloneCol, lwd = 3)
        #    invisible()
        #  })
        #}
      }
      
      if (cytoBand==TRUE){
        require(quantsmooth)
        par(xpd = NA)
        #paintCytobands(chrom=chr, units="bases", pos=c(0,(yrange[1]-0.5)), width=0.75, legend=F)	
      }
      
      if (!is.null(geneAnnot)){
        #par(xpd=F)
        colnames(geneAnnot) <- c("Gene","Chr","Start","Stop")
        geneAnnot <- geneAnnot[geneAnnot[,"Chr"]==chr,]
        if (nrow(geneAnnot) > 0){
        for (g in 1:dim(geneAnnot)[1]){
          print(geneAnnot[g,"Gene"])
          abline(v=as.numeric(geneAnnot[g,"Start"]),col="black",lty=3,xpd=F)
          abline(v=as.numeric(geneAnnot[g,"Stop"]),col="black",lty=3,xpd=F)			
          atP <- (as.numeric(geneAnnot[g,"Stop"]) - as.numeric(geneAnnot[g,"Start"]))/2 + as.numeric(geneAnnot[g,"Start"])
          if (atP < dataByChr[1,"Start"]){ atP <- dataByChr[1,"Start"] }
          else if (atP > dataByChr[dim(dataByChr)[1],"Start"]){ atP <- dataByChr[dim(dataByChr)[1],"Start"] }
          mtext(geneAnnot[g,"Gene"],side=3,line=0,at=atP,cex=gene.cex)
          }
        }
      }
    }
  }
}


plotIchorCNA <- function(dataIn, param = NULL, colName = "copy", segs=NULL, chr=NULL, ploidy = NULL, geneAnnot=NULL, yrange=c(-4,6), yaxis = "logRatio", xlim=NULL, xaxt = "n", cex = 0.5, gene.cex = 0.5, plot.title = NULL, cnCol = NULL, spacing=4, cytoBand=T, alphaVal=1, main){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  alphaSubcloneVal <- ceiling(alphaVal / 2 * 255); class(alphaVal) = "hexmode"
  subcloneCol <- "black" #c("#00FF00")
  if (is.null(cnCol)){
		cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 26))
	}
	cnCol <- paste(cnCol,alphaVal,sep="")
  names(cnCol) <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
#  segCol <- cnCol
#  ## add in colors for subclone if param provided
#  if (!is.null(param)){
#    ind <- ((which.max(param$ct) + 1) : length(param$ct)) + 1
#    cnCol[ind] <- paste0(cnCol[ind], alphaSubcloneVal / 2)
#    segCol[ind] <- "#00FF00"
#  }
  # adjust for ploidy #
  if (!is.null(ploidyT) & yaxis != "integer"){
    ploidyS <- purity * ploidy + (1-purity) * normCN
    dataIn[, colName] <- as.numeric(dataIn[, colName]) + log2(ploidy / 2)
    
    if (!is.null(segs) & yaxis != "integer"){
      segs[, colName] <- segs[, colName] + log2(ploidy / 2)
    }
  }
  dataIn <- dataIn[!is.na(dataIn[, colName]), ]
  
  if (!is.null(chr)){
    for (i in chr){
      
      dataByChr <- dataIn[dataIn[,"chr"]==as.character(i),]
          
       ## set y axis labels as either integer or logR copy number
       zero <- 0.25  
      cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      #ploidyToUse <- ploidyS
      if (i == "X"){
        normCN <- 1
        zero <- 0.5
        cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      }      
      if (yaxis == "integer"){
        y.ticks <- log2(cn)
        y.ticks[1] <- log2(zero)  
        yrange[1] <- y.ticks[1]    
        ylab <- "Copy Number"
        #dataByChr[, colName] <- log2(logRbasedCN(dataByChr[, colName], purity, ploidyT, cn=normCN))
        dataByChr[, colName] <- log2(dataByChr[, colName])
        if (!is.null(segs)){
      		segs[, colName] <- log2(segs[, colName])# + log2(ploidyS / 2)
    	}
        centreLine <- log2(normCN)
      }else{      
      	#dataByChr[, colName] <- dataByChr[, colName] + log2(ploidyS / 2)
      	cnLog <- log2(cn[-which(cn==3)] / normCN)  
        cn <- seq(-2,yrange[2],2)#c(-2, cn)
        y.ticks <- cn
        ylab <- "Copy Number (log2 ratio)"
        centreLine <- 0
      }
      
      # normCN <- round(ploidy)
#       cn <- c(0, 1:4, `^`(2, 3:(yrange[2]+1)))
#       ploidyToUse <- ploidyS
#       if (i == "X"){
#         normCN <- 1
#         ploidyToUse <- ploidyS / 2
#       }
#       if (yaxis == "integer"){
#         cnCor <- (purity*cn + (1-purity) * normCN)
#         #cnCor <- (normCN*purity*cn + normCN*(1-purity)*2) / ploidyToUse
#         cnCor.logR <- log2(cnCor/normCN)
#         y.ticks <- cnCor.logR
#         y.ticks[1] <- -2
#         ylab <- "Copy Number"
#       }else{
#         cn <- log2(cn[-which(cn==3)] / normCN)
#         cn[1] <- -2
#         y.ticks <- cn
#         cn <- format(cn, digits=2) 
#         ylab <- "Copy Number (log2 ratio)"
#       }

      #plot the data
      #if (outfile!=""){ pdf(outfile,width=10,height=6) }
      par(mar=c(spacing,8,4,2))
      #par(xpd=NA)
      coord <- (as.numeric(dataByChr[,"end"]) + as.numeric(dataByChr[,"start"]))/2
      if (is.null(xlim)){
        xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"start"]))
        xaxt <- "n"
      }
      if (is.null(plot.title)){
        plot.title <- paste("Chromosome ",i,sep="")
      }
      ## plot logR for bins ##
      plot(coord,as.numeric(dataByChr[, colName]),col=cnCol[dataByChr[,"event"]],
           pch=16, ylim=yrange, yaxt="n", 
           xlim=xlim, xaxt = xaxt, xlab="",ylab=ylab,
           cex.lab=1.5,cex.axis=1.5, cex=cex,las=1)
      axis(2, at=y.ticks, labels=cn, las=2, cex.axis=1.5)
      title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)
      ## plot centre line ##
      lines(c(1,as.numeric(dataByChr[dim(dataByChr)[1],3])),rep(0,2),type="l",col="grey",lwd=0.75)
      if (!is.null(segs)){
        segsByChr <- segs[segs[,"chr"]==as.character(i),,drop=FALSE]
        ind <- segsByChr$subclone.status == FALSE
        apply(segsByChr[ind, ], 1, function(x){
          lines(x[c("start","end")], rep(x[colName], 2), col = cnCol[x["call"]], lwd = 3)
          invisible()
        })
        if (sum(!ind) > 0){
          apply(segsByChr[!ind, ], 1, function(x){
            lines(x[c("start","end")], rep(x[colName], 2), col = subcloneCol, lwd = 3)
            invisible()
          })
        }
      }
      
      if (cytoBand==TRUE){
        require(quantsmooth)
        par(xpd = NA)
        #paintCytobands(chrom=chr, units="bases", pos=c(0,(yrange[1]-0.5)), width=0.75, legend=F)	
      }
      
      if (!is.null(geneAnnot)){
        #par(xpd=F)
        colnames(geneAnnot) <- c("Gene","Chr","Start","Stop")
        geneAnnot <- geneAnnot[geneAnnot[,"Chr"]==chr,]
        if (nrow(geneAnnot) > 0){
        for (g in 1:dim(geneAnnot)[1]){
          print(geneAnnot[g,"Gene"])
          abline(v=as.numeric(geneAnnot[g,"Start"]),col="black",lty=3,xpd=F)
          abline(v=as.numeric(geneAnnot[g,"Stop"]),col="black",lty=3,xpd=F)			
          atP <- (as.numeric(geneAnnot[g,"Stop"]) - as.numeric(geneAnnot[g,"Start"]))/2 + as.numeric(geneAnnot[g,"Start"])
          if (atP < dataByChr[1,"start"]){ atP <- dataByChr[1,"start"] }
          else if (atP > dataByChr[dim(dataByChr)[1],"start"]){ atP <- dataByChr[dim(dataByChr)[1],"start"] }
          mtext(geneAnnot[g,"Gene"],side=3,line=0,at=atP,cex=gene.cex)
          }
        }
      }
    }
  }else{  #plot for all chromosomes
    par(mar=c(spacing,8,2,2))
    midpt <- (as.numeric(dataIn[,"end"]) + as.numeric(dataIn[,"start"]))/2
    coord <- getGenomeWidePositions(dataIn[,"chr"],midpt)
    plot(coord$posns,as.numeric(dataIn[, colName]),
         col=cnCol[as.character(dataIn[,"event"])],pch=16,xaxt="n", ylim=yrange,
         xlim=c(1,as.numeric(coord$posns[length(coord$posns)])),
         xlab="",ylab="Copy Number (log2 ratio)",
         cex.lab=1.5,cex.axis=1.5,cex=0.5,las=1,bty="n",
         #main=dataIn[1,"sample"])
         main=main)
    #plot segments
    if (!is.null(segs)){
      coordEnd <- getGenomeWidePositions(segs[, "chr"], segs[, "end"])
      coordStart <- coordEnd$posns - (segs[, "end"] - segs[, "start"] + 1)
      xlim <- as.numeric(c(1, coordEnd$posns[length(coordEnd$posns)]))
      col <- cnCol[as.numeric(segs[, "state"] + 1)]
      value <- as.numeric(segs[, "median"])
      sc.status <- as.logical(segs[, "subclone.status"])
      mat <- as.data.frame(cbind(coordStart, coordEnd$posns, value, sc.status, col))
      rownames(mat) <- 1:nrow(mat)
      ## clonal CN
      ind <- mat$sc.status == FALSE
      apply(mat[ind, ], 1, function(x){
        lines(x[1:2], rep(x[3], 2), col = x[5], lwd = 3)
        invisible()
      })
      ## subclonal CN
      if (sum(!ind) > 0){
        apply(mat[!ind, ], 1, function(x){
          lines(x[1:2], rep(x[3], 2), col = subcloneCol, lwd = 3)
          invisible()
        })
      }
    }
    lines(as.numeric(c(1,coord$posns[length(coord$posns)])),rep(0,2),type="l",col="grey",lwd=2)
    plotChrLines(dataIn[,"chr"],coordEnd$chrBkpt,yrange)
  }
}

plotCNlogRByChr <- function(dataIn, colName = "copy", segs=NULL, chr=NULL, ploidy = NULL, geneAnnot=NULL, yrange=c(-4,6), xlim=NULL, cnCol = NULL, xaxt = "n", yaxis = "log2", cex = 0.5, gene.cex = 0.5, plot.title = NULL, spacing=4, cytoBand=T, alphaVal=1){
	#color coding
	alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
	if (is.null(cnCol)){
		cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 16))
		#cnCol <- col2rgb(c("green","darkgreen","blue","darkred","red","brightred"))
		names(cnCol) <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:15))
	}
	cnCol <- paste(cnCol,alphaVal,sep="")
	# adjust for ploidy #
	if (!is.null(ploidy)){
    dataIn[, colName] <- as.numeric(dataIn[, colName]) + log2(ploidy / 2)
    
    if (!is.null(segs)){
			segs[, "median"] <- segs[, "median"] + log2(ploidy / 2)
		}
  }
  

	if (!is.null(chr)){
	for (i in chr){
	
		dataByChr <- dataIn[dataIn[,"chr"]==as.character(i),]
		
		#plot the data
		#if (outfile!=""){ pdf(outfile,width=10,height=6) }
		par(mar=c(spacing,8,4,2))
		#par(xpd=NA)
		coord <- (as.numeric(dataByChr[,"end"]) + as.numeric(dataByChr[,"start"]))/2
		if (is.null(xlim)){
			xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"start"]))
			xaxt <- "n"
		}
		cn <- c(0, 1:4, `^`(2, 3:(yrange[2]+1)))
    normCN <- 2
    if (i == "X"){
    	normCN <- 1
    }
		if (yaxis == "integer"){
			y.ticks <- log2(cn / normCN)
			y.ticks[1] <- -2
			ylab <- "Copy Number"
		}else{
		  cn <- log2(cn[-which(cn==3)] / normCN)
		  cn[1] <- -2
			y.ticks <- cn
			cn <- format(cn, digits=2) 
			ylab <- "Copy Number (log2 ratio)"
		}
		if (is.null(plot.title)){
			plot.title <- paste("Chromosome ",i,sep="")
		}
    ## plot logR for bins ##
		plot(coord,as.numeric(dataByChr[, colName]),col=cnCol[as.numeric(dataByChr[,"state"])],
				pch=19, ylim=yrange,
				xlim=xlim, xaxt = xaxt, xlab="", yaxt="n",ylab=ylab,
				cex.lab=1.5,cex.axis=1.5, cex=cex,las=1)
		axis(2, at=y.ticks, labels=cn, las=2, cex.axis=1.5)
		title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)
    ## plot centre line ##
		lines(c(1,as.numeric(dataByChr[dim(dataByChr)[1],3])),rep(0,2),type="l",col="grey",lwd=0.75)
		if (!is.null(segs)){
			segsByChr <- segs[segs[,"chr"]==as.character(i),,drop=FALSE]
			tmp <- apply(segsByChr, 1, function(x){
			  lines(x[c("start","end")], rep(x["median"], 2), col = cnCol[as.numeric(x["state"])], lwd = 3)
			})
		}
    
		if (cytoBand==TRUE){
			par(xpd = NA)
			paintCytobands(chrom=chr, units="bases", pos=c(0,(yrange[1]-0.75)), width=0.5, legend=F)	
		}

		if (!is.null(geneAnnot)){
			#par(xpd=F)
			colnames(geneAnnot) <- c("Gene","Chr","Start","Stop")
			geneAnnot <- geneAnnot[geneAnnot[,"Chr"]==chr,]
			if (nrow(geneAnnot) > 0){
        for (g in 1:dim(geneAnnot)[1]){
          print(geneAnnot[g,"Gene"])
          abline(v=as.numeric(geneAnnot[g,"Start"]),col="black",lty=3,xpd=F)
          abline(v=as.numeric(geneAnnot[g,"Stop"]),col="black",lty=3,xpd=F)			
          atP <- (as.numeric(geneAnnot[g,"Stop"]) - as.numeric(geneAnnot[g,"Start"]))/2 + as.numeric(geneAnnot[g,"Start"])
          if (atP < dataByChr[1,"start"]){ atP <- dataByChr[1,"start"] }
          else if (atP > dataByChr[dim(dataByChr)[1],"start"]){ atP <- dataByChr[dim(dataByChr)[1],"start"] }
          mtext(geneAnnot[g,"Gene"],side=3,line=0,at=atP,cex=gene.cex)
        
        }
      }
		}
		}	
	  }else{  #plot for all chromosomes
	  	par(mar=c(spacing,8,2,2))
	  	midpt <- (as.numeric(dataIn[,"end"]) + as.numeric(dataIn[,"start"]))/2
    	coord <- getGenomeWidePositions(dataIn[,"chr"],midpt)
    	plot(coord$posns,as.numeric(dataIn[,colName]),
    			col=cnCol[as.numeric(dataIn[,"state"])],pch=19,xaxt="n", ylim=yrange,
    			xlim=c(1,as.numeric(coord$posns[length(coord$posns)])),
    			xlab="",ylab="Copy Number (log2 ratio)",
    			cex.lab=1.5,cex.axis=1.5,cex=0.5,las=1,bty="n",
    			main=dataIn[1,"sample"])
    	#plot segments
    	if (!is.null(segs)){
				coordEnd <- getGenomeWidePositions(segs[, "chr"], segs[, "end"])
				coordStart <- coordEnd$posns - (segs[, "end"] - segs[, "start"] + 1)
				xlim <- as.numeric(c(1, coordEnd$posns[length(coordEnd$posns)]))
				col <- cnCol[as.numeric(segs[, "state"])]
				value <- as.numeric(segs[, "median"])
				mat <- as.data.frame(cbind(coordStart, coordEnd$posns, value, col))
				rownames(mat) <- 1:nrow(mat)
				tmp <- apply(mat, 1, function(x){
					lines(x[1:2], rep(x[3], 2), col = x[4], lwd = 3)
				})
    	}
    	lines(as.numeric(c(1,coord$posns[length(coord$posns)])),rep(0,2),type="l",col="grey",lwd=2)
    	plotChrLines(dataIn[,"chr"],coord$chrBkpt,yrange)
    }
}

## modified to work for combine_TITAN_ICHOR/titan_ichor_cn.txt
findNearestLogR <- function(x, y, buffer = 1e6){
	y1 <- 0; y2 <- 0
	#y <- na.omit(y)
	bin1Ind <- x[["chromosome_1"]] == y[, "Chr"] & 
							x[["start_1"]] >= (y[, "Start"]) & 
							x[["start_1"]] <= (y[, "End"] + buffer)
	bin2Ind <- x[["chromosome_2"]] == y[, "Chr"] & 
							x[["start_2"]] >= (y[, "Start"]) & 
							x[["start_2"]] <= (y[, "End"] + buffer)
	if (sum(bin1Ind, na.rm=T) > 0){
	  # look at closest left and right points, find minimum distance to points, assign copy of min point
	  #ind <- (tail(which(bin1Ind), 1)-1):(tail(which(bin1Ind), 1))
	  ind <- tail(which(bin1Ind), 1)
	  #indClose <- ind[which.min(c(abs(y[ind[1], "end"]-x[["start_1"]]), abs(y[ind[2], "start"]-x[["start_1"]])))]
	  #y1 <- y[indClose, "copy"]
		#y1 <- max(y[tail(which(bin1Ind), 1):(tail(which(bin1Ind), 1)+1), "copy"], na.rm = TRUE)
		if (x[["orient_1"]] == "rev"){ # region to left of breakpoint 
      y1 <- y[max(ind - 1, 1), "LogRatio"]
    }else if (x[["orient_1"]] == "fwd"){  # region to right of breakpoint 
      y1 <- y[min(ind + 1, length(bin1Ind)), "LogRatio"]
    }
	}
	if (sum(bin2Ind, na.rm=T) > 0){
	  #ind <- (tail(which(bin2Ind), 1)-1):(tail(which(bin2Ind), 1))
	  ind <- tail(which(bin2Ind), 1)
	  #indClose <- ind[which.min(c(abs(y[ind[1], "end"]-x[["start_1"]]), abs(y[ind[2], "start"]-x[["start_1"]])))]
	  #y2 <- y[indClose, "copy"]
		#y2 <- max(y[tail(which(bin2Ind), 1):(tail(which(bin2Ind), 1)+1), "copy"], na.rm = TRUE)
		if (x[["orient_2"]] == "rev"){  # region to left of breakpoint 
      y2 <- y[max(ind - 1, 1), "LogRatio"]
    }else if (x[["orient_2"]] == "fwd"){  # region to right of breakpoint 
      y2 <- y[min(ind + 1, length(bin1Ind)), "LogRatio"]
    }
	}
	return(c(y1, y2))
}




curvedarrow.new <- function (from, to, lwd = 2, lty = 1, lcol = "black", arr.col = lcol,
    arr.pos = 0.5, curve = 1, dr = 0.01, endhead = FALSE, segment = c(0, 1), orient = "top", ...)
{
    dpos <- to - from
    angle <- atan(dpos[2]/dpos[1]) * 180/pi
    if (is.nan(angle))
        return
    mid <- 0.5 * (to + from)
    dst <- dist(rbind(to, from))
    ry <- abs(curve) #* dst
    aFrom <- 0
    aTo <- pi
    if (orient == "top") {
      aFrom <- 2 * pi
      aTo <- pi
    }else if (orient == "bottom"){
    	aFrom <- pi
    	aTo <- 2 * pi
    }
    if (segment[1] != 0)
        From <- segment[1] * aTo + (1 - segment[1]) * aFrom
    else From <- aFrom
    if (segment[2] != 1)
        To <- segment[2] * aTo + (1 - segment[2]) * aFrom
    else To <- aTo
    if ((from[1] < to[1] && orient == "top") || (from[1] > to[1] && orient == "bottom")){ ## draw arrow from down to upstream
      meanpi <- arr.pos * aFrom + (1 - arr.pos) * aTo
    }else{# if (from[1] > to[1]){ ## draw arrow from upstream to down 
      meanpi <- arr.pos * aTo + (1 - arr.pos) * aFrom
    }
    #if (endhead)
        #To <- meanpi 
    plotellipse(rx = dst/2, ry = ry, mid = mid, angle = angle,
        from = From, to = To, lwd = lwd, lty = lty, lcol = lcol)
    ell <- getellipse(rx = dst/2, ry = ry, mid = mid, angle = angle,
        from = 1.001 * meanpi, to = 0.999 * meanpi, dr = 0.002)
    if (endhead){
      if ((from[1] < to[1] && orient == "top") || (from[1] > to[1] && orient == "bottom")){ ## draw arrow from down to upstream
        code <- 2  # arrow draw at x1, y1
      }else{
        code <- 1 #arrow drawn at x0, y0
      }
      Arrows(ell[1, 1], y0=ell[1, 2], ell[nrow(ell), 1], y1=ell[nrow(ell), 2], 
          code = code, lcol = lcol, arr.col = arr.col, ...)
    }
    curvedarrow <- c(ell[nrow(ell), 1], ell[nrow(ell), 2])
}


## x is a GRangesList ##
# get x from load_snowman2
grl2circos <- function(x, region=NULL, genes = GRanges(), chr.exclude=NULL, cov=NULL, all.genes=FALSE, tracks.inside = 1, tracks.outside = 0, by.chromosome) {

  require(RCircos)
  if (is.null(region)) {
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
    #cyto.info <- data.frame(Chromosome=paste0("",names(seqlengths(si))), ChromStart=0, ChromEnd=as.numeric(seqlengths(si)), Band="p11", Stain="gvar", PlotColor="black")
    cyto.info <- cyto.info[!cyto.info$Chromosome %in% c("M","Y"),]
    if (!is.null(chr.exclude))
      cyto.info <- cyto.info[!cyto.info$Chromosome %in% chr.exclude, ]
  } else {

    rr <- .scale.gr(region, region)
    cyto.info <- data.frame(Chromosome=as.character(seqnames(rr)),
                            ChromStart=0,
                            ChromEnd=end(rr),
                            Band="p11",
                            Stain="gvar")

    ## restrict to just this region
    grl.in.region <- x[grl.in(x, region, only=TRUE)]
    if (length(grl.in.region)) {
      tmp <- .scale.gr(grl.unlist(grl.in.region), region)
      x <- split(tmp, tmp$uid)
    } else {
      x <- GRangesList()
    }
  }

  RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);
  rcircos.params <- RCircos.Get.Plot.Parameters();
  rcircos.params$base.per.unit <- 1e3
  rcircos.params$track.height <- 0.2
  RCircos.Reset.Plot.Parameters(rcircos.params)

  RCircos.Set.Plot.Area()

  #track.num=1
  #rr <- gr2dt(region)
  #setnames(rr, c("seqnames","start","end"), c("Chromosome","chromStart","chromEnd"))
  #RCircos.Gene.Connector.Plot(rr, track.num, "out");
  #rr$names <- LETTERS[1:4]
  #track.num=2
  #name.col=5
  #RCircos.Gene.Name.Plot(rr, name.col,track.num,"out");
  #RCircos.Label.Chromosome.Names()
  
  ## get the gene label dat
  gene.dat <- data.frame()
  if (length(genes)) {
    genes <- .scale.gr(genes,region)
    ## restrict to just this region

    if (!all.genes) {
      gr1 <- grl.unlist(x)[grl.unlist(x)$grl.iix==1]
      gr2 <- grl.unlist(x)[grl.unlist(x)$grl.iix==2]    
      fo1 <- gr.findoverlaps(gr1+10e3, genes)
      fo2 <- gr.findoverlaps(gr2+10e3, genes)
      
      fo <- c(fo1,fo2)
      gene.dat <- data.frame()
      if (length(fo)) {
        fo <- fo[!duplicated(fo$subject.id)]
        gene.dat = data.frame(Chromosome = seqnames(genes[fo$subject.id]), chromStart=start(genes[fo$subject.id]),
          chromEnd=end(genes[fo$subject.id]), Gene=genes$gene[fo$subject.id])
      }
    } else {
      gene.dat = data.frame(Chromosome = seqnames(genes), chromStart=start(genes),
        chromEnd=end(genes), Gene=genes$gene)
    }

    
  }
  gename.col <- 4;

  ## set the links data
  links <- grl2links(x)

  ## draw
  RCircos.Chromosome.Ideogram.Plot();
  if (nrow(gene.dat) > 0) {
    track.num=2
    RCircos.Gene.Connector.Plot(gene.dat, track.num, side);
    track.num <- 1;
    name.col <- 4;
    RCircos.Gene.Name.Plot(gene.dat, name.col,track.num, side);
    track.num <- 3
    RCircos.Tile.Plot(gene.dat, track.num, "in");
  }

  ## coverage plot if there
  if (!is.null(cov)) {
    cov <- .scale.gr(cov, region)
    ##covd <- gr.val(gr.tile(reduce(cov), w=), cov, val="ratio")
    track.num=4
    by.fold = 1
    ttt <- gr2dt(cov)
    setnames(ttt, c("seqnames","end"), c("chromosome","stop"))
    data.col <- match("ratio.normalized", colnames(ttt))
    RCircos.Scatter.Plot(as.data.frame(ttt[ratio.normalized > 0 & ratio.normalized < 4]), data.col,
                         track.num, "in", by.fold=by.fold, min.value=0, max.value=3)
  }
  
  side <- "in";
  track.num <- 1;
  if (nrow(links) > 0)
  	#RCircos.Get.Link.Colors()
    RCircos.Link.Plot(links, track.num, by.chromosome=by.chromosome) ## by.chromosome is for color

}

load_snowman2 <- function(x) {
  wc <- as.numeric(system(paste("grep -v ^#", x, "| wc -l"), intern=TRUE))
  if (!wc)
    return (GRangesList());
  ff <- fread(paste("grep -v ^#", x), sep="\t")
  setnames(ff, c("V1","V2","V3","V4","V5","V6","V7", "V8"), c("seqnames", "start","id","REF", "ALT", "QUALITY","FILTER","INFO"))
  ff$sample <- basename(x)
  ff[, uid := paste(id, sample, sep="_")]
  ff[, per_sample_count := length(uid)/2, by=sample]
  ff[, ID := gsub("([0-9]+):(1|2)", "\\1", id)]
  ff[, EVDNC := gsub(".*?EVDNC=([A-Z_]+).*", "\\1", INFO)]
  ff[, NUMPARTS := as.integer(gsub(".*?NUMPARTS=([0-9]+).*", "\\1", INFO))]
  ff[, SOMATIC := grepl("SOMATIC", INFO)]
  ff[, SCTG := gsub(".*?SCTG=(.*?);.*", "\\1", INFO)]
  ff[, DISC_MAPQ := as.numeric(gsub(".*?DISC_MAPQ=([0-9]+).*", "\\1", INFO))]
  ff[, NM := as.integer(gsub(".*?;NM=([0-9]+).*", "\\1", INFO))]
  ff[, MATENM := as.integer(gsub(".*?;MATENM=([0-9]+).*", "\\1", INFO))]
  ff[, MAPQ := as.integer(gsub(".*?;MAPQ=([0-9]+).*", "\\1", INFO))]
  print("...still formatting")
  ff[, REPSEQ := gsub(".*?;REPSEQ=([A-Z]+).*", "\\1", INFO)]
  ff[, REPSEQ := ifelse(grepl(";", REPSEQ), "", REPSEQ)] 
  ff[, HOMSEQ := gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO)] ##{ xx <- gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO); if (!grepl(";",xx)) { xx } else { "" }}]
  ff[, HOMSEQ := ifelse(grepl(";", HOMSEQ), "", HOMSEQ)] 
  ff[, INSERTION := gsub(".*?;INSERTION=([A-Z]+).*", "\\1", INFO)] ##{ xx <- gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO); if (!grepl(";",xx)) { xx } else { "" }}]
  ff[, INSERTION := ifelse(grepl(";", INSERTION), "", INSERTION)]

  if (c("V11") %in% colnames(ff)) {
    ff$TUMALT = sapply(strsplit(ff$V11, ":"), function(y) {as.numeric(y[2])})
    ff$TUMLOD = sapply(strsplit(ff$V11, ":"), function(y) {as.numeric(y[9])})
  }
  if ("V10" %in% colnames(ff)) {
    ff$NORMALT = sapply(strsplit(ff$V10, ":"), function(y) {as.numeric(y[2])})
    ff$NORMLOD = sapply(strsplit(ff$V10, ":"), function(y) {as.numeric(y[9])})
  }
  ff[, strand := ifelse(grepl("^\\[", ALT) | grepl("^\\]", ALT), '-', '+')]
  ff[, inv := strand[1] == strand[2], by=c("sample","ID")]
  ff[, altstrand := rev(strand), by=c("sample", "ID")]
  ff[, altpos := as.integer(gsub(".*?:([0-9]+).*", "\\1", ALT))]
  ff[, altchr := gsub(".*?(\\[|\\])(.*?):([0-9]+).*", "\\2", ALT)]
  ff[, SPAN := ifelse(seqnames==altchr, abs(start - altpos), -1)]
  
  ff[, c("INFO","id","V9", "V10") := NULL]
  setkey(ff, ID)
  gr <- with(ff, GRanges(seqnames, IRanges(start, start), strand=strand, ID=ID))
  grl <- split(gr, gr$ID)
  fff <- as.data.frame(ff)
  fff <- fff[!duplicated(fff$ID), setdiff(colnames(fff), c("seqnames","start","strand"))]
  mcols(grl) <- fff
  return(grl)
}

grl2links <- function(x) {

  y = grl.unlist(x);

  ix = y$grl.iix == 1;
  sn = as.character(seqnames(y))
  p = as.numeric(start(y));
  links = data.frame(Chromosome=as.character(sn[ix]), chromStart=p[ix], chromEnd=p[ix], Chromsome.1=as.character(sn[!ix]), chromStart.1=p[!ix], chromEnd.1=p[!ix], stringsAsFactors=FALSE)
  return (links)
  
}


