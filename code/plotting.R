###########################################
############ PLOTTING FUNCTION ###########
###########################################
## sv is a data.table output by combineSVABAandTITAN.R
### TO-DO: 
# 1) Filter by support
# 2) Strand direction of breakpoint
plotRearrangementArcs <- function (sv, cn, ploidy = NULL, interchr=TRUE, xlim=NULL, segIn=NULL, support=1, chr=NULL, chrLens=NULL, minSPAN = 10, buffer=1e6, plotAtCentre = FALSE, centreLine=0, arcHeight=4, lty=1, lcol = "black", svTypeCol=FALSE, arr.col = "black", arr.pos = 1, endhead = FALSE, lwd = 1, orient="topbottom", include.inter.chr = FALSE, offset.factor = 1.15){
	require(diagram)
	sv <- sv[sv$SPAN >= minSPAN | sv$SPAN == -1, ]
	sv <- as.data.frame(sv)
	
	# adjust for ploidy #
	if (!is.null(ploidy)){
    cn[, "LogRatio"] <- as.numeric(cn[, "LogRatio"]) + log2(ploidy / 2)
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

		## plot intra-chromosomal SV ##
		
			for (i in 1:nrow(dataByChr)){
				if (dataByChr[i, "inWindow"] == "NONE"){
					next
				}
				#print(dataByChr[i, ])
				## get cn logR of nearest bin
				yLogR <- findNearestLogR(dataByChr[i, ], cn, buffer = buffer * 100)	
				yToUse <- yLogR			
				orient <- "top"
				start1 <- dataByChr[i,"start_1"]
        end1 <- dataByChr[i,"start_2"]

				if (plotAtCentre && !is.null(centreLine)){
					yToUse <- c(centreLine, centreLine)
				}
				curveHeight <- arcHeight[2] - mean(yToUse)#/ abs(start1 - end1)
				text.height <- arcHeight[2]
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
						text.height <- 1/arcHeight[2]			
						if (!is.null(centreLine)) { yToUse <- c(centreLine, centreLine) }
						curveHeight <- arcHeight[1] + mean(yToUse)
					}
	
				 	if (is.null(offset.factor)) { offset.factor <- 1.15 }
					leftBreak <- paste0(dataByChr[i,"chromosome_1"], ":", dataByChr[i, "start_1"])
					rightBreak <- paste0(dataByChr[i,"chromosome_2"], ":", dataByChr[i, "start_2"])
				 	## for inter chromosome
				 	interchr.offset <- (par("usr")[2] - par("usr")[1]) * offset.factor # left coordinate if outside plot
					if (dataByChr[i, "chromosome_1"] != chr){					    
						if (abs(dataByChr[i, "start_2"]-xlim[1]) < abs(dataByChr[i, "start_2"]-xlim[2])){
						# 2nd bkpt is closer to chr begin
							start1 <- par("usr")[1] - interchr.offset # left boundary  
							#text(x=par("usr")[1], y=text.height, pos=4, cex=0.5, label=leftBreak) 
						}else{
							start1 <- end1
							end1 <- par("usr")[2] + interchr.offset#right boundary  
							#text(x=par("usr")[2], y=text.height, pos=2, cex=0.5, label=rightBreak)
							#yToUse <- yToUse[c(2,1)]              
						}                   
					}else if (dataByChr[i, "chromosome_2"] != chr){
						if (abs(dataByChr[i, "start_1"]-xlim[1]) < abs(dataByChr[i, "start_1"]-xlim[2])){
						# 1st bkpt is closer to chr begin
							end1 <- start1
							start1 <- par("usr")[1] - interchr.offset #- (par("usr")[2] * offset.factor - par("usr")[2]) # left boundary 
							#text(x=par("usr")[1], y=text.height, pos=4, cex=0.5, label=leftBreak) 
							#yToUse <- yToUse[c(2,1)]               
						}else{
							end1 <- par("usr")[2] + interchr.offset #right boundary 
							#text(x=par("usr")[2], y=text.height, pos=2, cex=0.5, label=rightBreak)     
						}  
					}
					## for zoom within same chromosome
					if (dataByChr[i, "intraChr"]){
						if (start1 < xlim[1]){ 
							start1 <- par("usr")[1] - interchr.offset #- (par("usr")[2] * offset.factor - par("usr")[2]) # left boundary 
							#text(x=par("usr")[1], y=text.height, pos=4, cex=0.5, label=leftBreak) 
						}
						if (start1 > xlim[2]){
							start1 <- par("usr")[2] - interchr.offset #- (par("usr")[2] * offset.factor - par("usr")[2]) #right boundary
						}
						if (end1 < xlim[1]){
							end1 <- par("usr")[1] + interchr.offset # left boundary
						}
						if (end1 > xlim[2]){
							end1 <- par("usr")[2] + interchr.offset # left boundary   
							#text(x=par("usr")[2], y=text.height, pos=2, cex=0.5, label=rightBreak)             
						}        
						#if (sv.type == "inv") { end2 <- start1; start2 <- end1 }    
					}
					#curveHeight <- par("usr")[3]# - mean(yToUse)#* 2 / (end1 - start1)
				}
				
				curvedarrow.new(from = c(start1, yToUse[1]), endhead = endhead, orient = orient,
									to = c(end1, yToUse[2]), lty = lty, lwd = lwd, lcol=lcol.toUse,
									arr.col=arr.col.toUse, arr.width = 0.1, arr.length=0.2,
									curve = curveHeight, arr.type = "triangle", arr.pos = arr.pos)

			}
		}		
	}
	invisible()
}

## format of dataIn is output from /home/unix/gavinha/software/code/git/scripts/titan/analysis/combineTITAN-ichor.R
plotTitanIchorCNA <- function(dataIn, param = NULL, colName = "LogRatio", callColName="Corrected_Call", segs=NULL, chr=NULL, purity = NULL, ploidyT = NULL, geneAnnot=NULL, yrange=c(-4,6), yaxis = "logRatio", xlim=NULL, xaxt = "n", cex = 0.5, gene.cex = 0.5, plot.title = NULL, cnCol = NULL, spacing=4, cytoBand=T, alphaVal=1, main){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  alphaSubcloneVal <- ceiling(alphaVal / 2 * 255); class(alphaVal) = "hexmode"
  subcloneCol <- c("#00FF00")
  if (is.null(cnCol)){
		cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 26))
		cnCol <- c(cnCol, "HET"="#0000FF", "DLOH"="#006400", "NLOH"="#0000FF", "ALOH"="#FF0000", "ASCNA"="#FF0000", "BCNA"="#FF0000", "UBCNA"="#FF0000")
	}else{
		cnCol.col <- as.character(cnCol[1])
		cnCol <- c(cnCol, "HET"=cnCol.col, "DLOH"=cnCol.col, "NLOH"=cnCol.col, "ALOH"=cnCol.col, "ASCNA"=cnCol.col, "BCNA"=cnCol.col, "UBCNA"=cnCol.col)
	}
	names(cnCol)[1:30] <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
	#cnCol <- paste(cnCol,alphaVal,sep="")
  # adjust for ploidy #
  normCN <- 2
  if (!is.null(ploidyT) & yaxis != "integer"){
    ploidyS <- purity * ploidyT + (1-purity) * normCN
    dataIn[, colName] <- as.numeric(dataIn[, colName]) + log2(ploidyS / 2)
    
    if (!is.null(segs)){
      segs[, colName] <- segs[, colName] + log2(ploidyS / 2)
    }
  }
  
  
  if (!is.null(chr)){
    for (i in chr){
      dataByChr <- dataIn[dataIn[,"Chr"]==as.character(i),]
       ## set y axis labels as either integer or logR copy number
      #avgTumPloidy <- round(ploidyT)
 
      zero <- 0.5  
      cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      #ploidyToUse <- ploidyS
      if (i == "X"){
        normCN <- 1
        zero <- 0.25
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
      plot(coord,as.numeric(dataByChr[, colName]),col=cnCol[dataByChr[,callColName]],
           pch=16, ylim=yrange, yaxt="n",
           xlim=xlim, xaxt = xaxt, xlab="",ylab=ylab,
           cex.lab=1.5,cex.axis=1.5, cex=cex,las=1)
      axis(2, at=y.ticks, labels=cn, las=2, cex.axis=1.5)
      title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)
      ## plot centre line ##
      lines(c(1,tail(na.omit(dataByChr[,3]), 1)),rep(centreLine,2),type="l",col="grey",lwd=0.75)
      if (!is.null(segs)){
        segsByChr <- segs[segs[,"Chromosome"]==as.character(i),,drop=FALSE]
        #ind <- segsByChr$subclone.status == FALSE
        apply(segsByChr, 1, function(x){
          lines(x[c("Start","End")], rep(x[colName], 2), col = cnCol[x[callColName]], lwd = 3)
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

## compute copy number using corrected log ratio ##
logRbasedCN <- function(x, purity, ploidyT, cn = 2){
	ct <- (2^x * (cn * (1 - purity) + purity * ploidyT * (cn / 2)) - cn * (1 - purity)) / purity
	ct <- sapply(ct, max, 1/2^6)
	return(ct)
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
      #y1 <- tail(na.omit(y[max((ind-10):ind), "LogRatio"]), 1)
    }else if (x[["orient_1"]] == "fwd"){  # region to right of breakpoint 
      y1 <- y[min(ind + 1, length(bin1Ind)), "LogRatio"]
      #y1 <- head(na.omit(y[max(ind:(ind+10)), "LogRatio"]), 1)
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
	if (is.na(y1)) { y1 <- 0 }
	if (is.na(y2)) { y2 <- 0 }
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

### modify SNPchip function "plotIdiogram"
plotIdiogram.hg38 <- function (chromosome, cytoband, seqinfo, cytoband.ycoords, xlim,
                               ylim = c(0, 2), new = TRUE, label.cytoband = TRUE, label.y = NULL,
                               srt, cex.axis = 1, outer = FALSE, taper = 0.15, verbose = FALSE,
                               unit = c("bp", "Mb"), is.lattice = FALSE, ...)
{
  def.par <- par(no.readonly = TRUE, mar = c(4.1, 0.1, 3.1,
                                             2.1))
  on.exit(def.par)
  if (is.lattice) {
    segments <- lsegments
    polygon <- lpolygon
  }
  
  cytoband <- cytoband[cytoband[, "chrom"] == chromosome, ]
  unit <- match.arg(unit)
  if (unit == "Mb") {
    cytoband$start <- cytoband$start/1e+06
    cytoband$end <- cytoband$end/1e+06
  }
  if (missing(cytoband.ycoords)) {
    cytoband.ycoords <- ylim
  }
  rownames(cytoband) <- as.character(cytoband[, "name"])
  sl <- seqlengths(seqinfo)[chromosome]
  if (missing(xlim))
    xlim <- c(0, sl)
  if (unit == "Mb")
    xlim <- xlim/1e+06
  cytoband_p <- cytoband[grep("^p", rownames(cytoband), value = TRUE),
                         ]
  cytoband_q <- cytoband[grep("^q", rownames(cytoband), value = TRUE),
                         ]
  p.bands <- nrow(cytoband_p)
  cut.left <- c()
  cut.right <- c()
  for (i in seq_len(nrow(cytoband))) {
    if (i == 1) {
      cut.left[i] <- TRUE
      cut.right[i] <- FALSE
    }
    else if (i == p.bands) {
      cut.left[i] <- FALSE
      cut.right[i] <- TRUE
    }
    else if (i == (p.bands + 1)) {
      cut.left[i] <- TRUE
      cut.right[i] <- FALSE
    }
    else if (i == nrow(cytoband)) {
      cut.left[i] <- FALSE
      cut.right[i] <- TRUE
    }
    else {
      cut.left[i] <- FALSE
      cut.right[i] <- FALSE
    }
  }
  for (i in seq_len(nrow(cytoband))) {
    if (as.character(cytoband[i, "gieStain"]) == "stalk") {
      cut.right[i - 1] <- TRUE
      cut.left[i] <- NA
      cut.right[i] <- NA
      cut.left[i + 1] <- TRUE
    }
  }
  include <- cytoband[, "end"] > xlim[1] & cytoband[, "start"] <
    xlim[2]
  cytoband <- cytoband[include, ]
  N <- nrow(cytoband)
  cut.left <- cut.left[include]
  cut.right <- cut.right[include]
  if (new) {
    xx <- c(0, cytoband[nrow(cytoband), "end"])
    yy <- cytoband.ycoords
    plot(xx, yy, xlim = xlim, type = "n", xlab = "", ylab = "",
         axes = FALSE, yaxs = "i", ylim = ylim, ...)
  }
  top <- cytoband.ycoords[2]
  bot <- cytoband.ycoords[1]
  h <- top - bot
  p <- taper
  for (i in seq_len(nrow(cytoband))) {
    start <- cytoband[i, "start"]
    last <- cytoband[i, "end"]
    delta = (last - start)/4
    getStain <- function(stain) {
      switch(stain, gneg = "grey100", gpos25 = "grey90",
             gpos50 = "grey70", gpos75 = "grey40", gpos100 = "grey0",
             gvar = "grey100", stalk = "brown3", acen = "brown4",
             "white")
    }
    color <- getStain(as.character(cytoband[i, "gieStain"]))
    if (is.na(cut.left[i]) & is.na(cut.right[i])) {
      delta <- (last - start)/3
      segments(start + delta, cytoband.ycoords[1], start +
                 delta, cytoband.ycoords[2])
      segments(last - delta, cytoband.ycoords[1], last -
                 delta, cytoband.ycoords[2])
    }
    else if (cut.left[i] & cut.right[i]) {
      yy <- c(bot + p * h, bot, bot, bot + p * h, top -
                p * h, top, top, top - p * h)
      polygon(c(start, start + delta, last - delta, last,
                last, last - delta, start + delta, start), yy,
              col = color)
    }
    else if (cut.left[i]) {
      yy <- c(bot + p * h, bot, bot, top, top, top - p *
                h)
      polygon(c(start, start + delta, last, last, start +
                  delta, start), yy, col = color)
    }
    else if (cut.right[i]) {
      yy <- c(bot, bot, bot + p * h, top - p * h, top,
              top)
      polygon(c(start, last - delta, last, last, last -
                  delta, start), yy, col = color)
    }
    else {
      polygon(c(start, last, last, start), c(bot, bot,
                                             top, top), col = color)
    }
  }
  my.x <- (cytoband[, "start"] + cytoband[, "end"])/2
  if (label.cytoband & !is.lattice) {
    if (is.null(label.y)) {
      axis(1, at = my.x, labels = rownames(cytoband), outer = outer,
           cex.axis = cex.axis, line = 1, las = 3, tick = FALSE)
      axis(1, at = cytoband$start, outer = outer, cex.axis = cex.axis,
           line = 1, las = 3, labels = FALSE)
    }
    else {
      if (!is.numeric(label.y)) {
        warning("label.y must be numeric -- using default y coordinates for cytoband labels")
        label.y <- bot - p * h
      }
      if (missing(srt))
        srt <- 90
      text(x = my.x, y = rep(label.y, length(my.x)), labels = rownames(cytoband),
           srt = srt)
    }
  }
  return()
}


