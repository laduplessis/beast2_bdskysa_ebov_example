library(bdskytools)
library(yaml)

dateorigin <- "1970-01-01"

#' TMRCA 95% HPD on current set of axes
plotTMRCA <- function(TreeHeight, xlim, ylim, col=pal.dark(cgreen,0.5)) { 
  mrcahpd     <- getHPD(TreeHeight)
  mrcadensity <- density(TreeHeight, from=mrcahpd[1], to=mrcahpd[3])
  scaling     <- 0.9*ylim[2]/max(mrcadensity$y)
  end         <- as.Date(xlim[2])
  #idxs        <- which(end-(365*mrcadensity$x) <= plotstart)
  mrcax       <- end-(365*mrcadensity$x)
  mrcay       <- mrcadensity$y*scaling
  
  usr <- par("usr")
  clip(xlim[1], xlim[2], ylim[1], ylim[2])
  
  polygon(c(max(mrcax), mrcax, min(mrcax)), c(0, mrcay, 0), col=col, border=NA)
  
  #text(0, -ReMax*0.018, paste0("Median TMRCA = ", format.Date(xmax-(mrca_hpd[2])*365, format="%Y/%m/%d")), 
  #    col=pal.dark(cgreen), pos=4, cex=0.8)
  #print(median(TreeHeight))
  #print(end-365*median(TreeHeight))  
  
  do.call("clip", as.list(usr))
}


#' Histogram of cases on current set of axes
plotCaseHist <- function(dataset, daterange, ylim=NA, ymax=NA, breaks='weeks', col=pal.dark(cblue), side=2, line=3, ylab="",
                         start.on.monday=FALSE, plotCounts=FALSE, lowCountMark=10, cex.label=1.2, cex.axis=0.8, cex.counts=0.6, ...) {
  
  caseslist <- c()
  i <- 1
  for (i in 1:nrow(dataset)) {
    caseslist <- c(caseslist, rep(rownames(dataset)[i], dataset$combined[i]))
    i <- i + 1
  }
  caseslist <- as.Date(caseslist)
  
  # Get weekly and monthly breaks across daterange
  # weeks  <- as.Date(hist(daterange, breaks='weeks',  plot=FALSE, plot.on.monday=start.on.monday)$breaks, origin=dateorigin)
  # months <- as.Date(hist(daterange, breaks='months', plot=FALSE)$breaks, origin=dateorigin)
  
  # Get counts in bins
  z <- hist(caseslist, breaks=breaks, freq=TRUE, plot = FALSE, start.on.monday=start.on.monday)
  
  # Get Y limits
  if (is.na(ymax)) {
    ymax <- max(z$counts)
  }
  maxprop   <- 0.8
  yticks    <- pretty(c(0,ymax/maxprop))
  ymax      <- max(yticks)
  scaling   <- ylim[2]/ymax
  abscounts <- z$counts
  z$counts  <- scaling*z$counts
  
  # Start the plot
  #plot(1,type='n',xlim=daterange, ylim=ylims, xaxs='i',yaxs='i',bty='n', xaxt='n', yaxt='n', ...)
  
  # Shade months
  #for (i in seq(1,length(months),by=2)) {
  #  rect(months[i],ylims[1], months[i+1], ylims[2], col = pal.light(cgray), border=NA)
  #}
  
  # Plot y-axis
  axis(side, at=(scaling*yticks), labels=yticks, las=2, cex.axis=cex.axis)
  mtext(side=side, line=line, text=ylab, cex=cex.label)
  
  
  # Don't plot x-axis
  #axis(1, at=months, labels=NA)
  #text(x=months[1:(length(months)-1)], y=-0.04*ylims[2], format.Date(months[2:length(months)],format="%b"), xpd=TRUE, srt=0, pos=4, cex=cex.months)
  #yearmids <- which(format.Date(months,format="%m") == "06")
  #text(x=months[yearmids], y=-0.06*ylims[2], format.Date(months[yearmids],"%Y"), xpd=TRUE, srt=0, pos=1, cex=cex.years)
  
  # Plot histogram bins
  lines(z, border=pal.dark(cwhite), col=col)
  
  # Add counts
  if (plotCounts) {
    countlabels <- abscounts
    print(countlabels)
    countlabels[which(countlabels==0)] <- NA
    text(z$mids, z$counts+(0.02*ylim[2]), labels=countlabels, col="red", cex=cex.counts, srt=90, adj=c(0,0.5),xpd=FALSE)
  }
  
  # Add arrows
  if (lowCountMark > 0) {
    lowCounts <- z$mids[which(abscounts < lowCountMark)]
    text(lowCounts, rep((scaling*lowCountMark)+(0.04*ylim[2]),length(lowCounts)), label=rep(expression(symbol("\257")), length(lowCounts)), col="red", cex=cex.axis)
  }
}



#' Plot reproductive number and sampling proportion skylines
plotSkyline <- function(configfile, logfile, samplingfile=NA, outfile, reported=NULL, 
                        burnin=0.1, 
                        startdate, enddate, Regrid=50, Remax=3, samplingmax=1,
                        plotMarks=c(), plotMRCA=TRUE, plotOrigin=TRUE) {
  
  # Setup
  pars <- yaml.load_file(configfile)
  lf   <- readLogfile(logfile, burnin=burnin)
  
  endyear   <- getYearDate(enddate)
  startyear <- getYearDate(startdate)
  
  
  # Sampling proportion
  if (grepl("treeslicer",logfile)) {
    samplingtimes <- c(startyear,unlist(rev(getSkylineSubset(lf, "SamplingTreeSlice.dates[0-9]")[1,]), use.names=FALSE))
  } else {
    samplingparts <- strsplit(pars$samplingRateChangeTimes,"\t")[[1]]
    if (length(samplingparts) == 1)
      samplingparts <- strsplit(pars$samplingRateChangeTimes," ")[[1]]  
    samplingtimes <- sort(c(startyear, endyear - as.numeric(samplingparts)))
  }
  samplingdates <- getDayDate(samplingtimes)
  sampling_sky  <- getSkylineSubset(lf, "samplingProportion.[0-9]")
  sampling_hpd  <- getMatrixHPD(sampling_sky, alpha=0.05)
  rownames(sampling_hpd) <- c("Lower","Median","Upper") 
  colnames(sampling_hpd) <- samplingdates[1:(length(samplingdates)-1)] 
  
  
  # Re
  Re_sky     <- getSkylineSubset(lf, "reproductiveNumber.[0-9]")
  if (grepl("treeslicer",logfile) || grepl("rootcondition",logfile)) {
    Re_gridded <- gridSkylineDates(Re_sky, lf$TreeHeight, enddate=enddate, from=startdate, intervals=Regrid)
  } else {
    Re_gridded <- gridSkylineDates(Re_sky, lf$origin, enddate=enddate, from=startdate, intervals=Regrid)
  }
  Re_gridded_hpd   <- getMatrixHPD(Re_gridded$skyline, alpha=0.05)
  colnames(Re_gridded_hpd) <- format.Date(Re_gridded$dates,origin=dateorigin)
  rownames(Re_gridded_hpd) <- c("Lower","Median","Upper")     
  
  
  # X-axis ticks
  xticks     <- getMonths(Re_gridded$dates[1],Re_gridded$dates[length(Re_gridded$dates)])
  start      <- findInterval(enddate-(365*median(lf$TreeHeight)), Re_gridded$dates)
  end        <- length(Re_gridded$dates)
  yearstarts <- which(format.Date(xticks,format="%m") == "01")
  monthlabs  <- format.Date(xticks,format="%b ")
  yearlabs   <- format.Date(xticks[yearstarts],"%Y")
  
  
  #NewFig(outfile,width=7.5, aspectratio = 5)
  #par(mar=c(3,5,2,4)+0.1, mfrow=c(1,2))
  
  
  ############################
  # Plot reproductive number #
  ############################
  NewFig(paste0(outfile,".ReproductiveNumber.pdf"),width=3.5, aspectratio = 5/2)
  par(mar=c(3,5,2,4)+0.1)
  
  #plotSkylinePretty(Re_gridded$dates[start:end], Re_gridded_hpd[,start:end], type='step',
  #                  col=NA, fill=NA, bty='o',
  #                  xaxis=FALSE, yaxis=FALSE, axispadding=0, 
  #                  xlims=range(xticks), ylims=c(0,Remax))
  plot(1, type='n',  xlim=range(xticks), ylim=c(0,Remax), xaxs='i', yaxs='i', xlab="", ylab="", axes=FALSE)
  
  # Shade months
  for (i in seq(1,length(xticks),by=2)) {
    rect(xticks[i],0, xticks[i+1],Remax, col = pal.light(cgray), border=NA)
  }
  
  # Add histogram of reported cases
  if (!is.null(reported)) {
      plotCaseHist(reported, range(xticks), ylim=c(0,Remax), 
                   side=4, ylab="Confirmed cases (weekly)", cex.label=1.0, cex.axis=0.7,
                   plotCounts=FALSE, lowCountMark=10)
  }
  
  # Add Re skyline
  plotSkylinePretty(Re_gridded$dates[start:end], Re_gridded_hpd[,start:end], type='smooth',
                    col=pal.dark(corange), fill=pal.dark(corange,0.7),
                    xaxis=TRUE, yaxis=TRUE, axispadding=0, ylims=c(0,Remax),
                    ylab=expression("R"[e]), xlab="", xline=2, side=2, yline=2.5,
                    lwd=1, cex.axis=0.7, cex.label=1.0, xticks=xticks, xticklabels=NA,
                    new=FALSE, add=TRUE)
  
  # Line at 1 and box around plot
  lines(range(xticks), rep(1,2), col=pal.dark(cred), lty=2, lwd=1)
  rect(min(xticks), 0, max(xticks), Remax, xpd=TRUE)
  
  # Add marks
  if (length(plotMarks) > 0) {
      for (i in 1:length(plotMarks)) {
        abline(v=plotMarks[i], lty=3)
        text(x=plotMarks[i], y=Remax, LETTERS[i], xpd=TRUE, srt=0, pos=3, cex=1.0)
      }
  }
    
  # Plot X-axis
  text(x=xticks[1:(length(xticks)-1)], y=-0.02*Remax, monthlabs[1:(length(monthlabs)-1)], xpd=TRUE, srt=0, pos=1, cex=0.7)
  text(x=xticks[yearstarts], y=-0.1*Remax, yearlabs, xpd=TRUE, srt=0, pos=1, cex=0.7)
  
  # Add origin/TMRCA
  if (plotMRCA == TRUE) plotTMRCA(lf$TreeHeight, range(Re_gridded$dates), c(0,Remax), col=pal.dark(cpurple, 0.7)) 
  
  if (plotOrigin == TRUE) {
    if ("origin" %in% names(lf))
      plotTMRCA(lf$origin, range(Re_gridded$dates), c(0,Remax), col=pal.dark(cgreen, 0.7)) 
  }
  
  dev.off()
  
  ############################
  # Plot sampling proportion #
  ############################
  NewFig(paste0(outfile, ".SamplingProportion.pdf"),width=3.5, aspectratio = 5/2)
  par(mar=c(3,5,2,4)+0.1)
  
  #plotSkylinePretty(samplingdates, sampling_hpd, type='step',
  #                  col=NA, fill=NA,
  #                  xaxis=FALSE, yaxis=FALSE, axispadding=0, 
  #                  xlims=range(xticks), ylims=c(0,samplingmax))
  plot(1, type='n',  xlim=range(xticks), ylim=c(0,samplingmax), xaxs='i', yaxs='i', xlab="", ylab="", axes=FALSE)
  
  # Shade months
  for (i in seq(1,length(xticks),by=2)) {
    rect(xticks[i],0, xticks[i+1],samplingmax, col = pal.light(cgray), border=NA)
  }
  
  # Add sampling proportion skyline
  plotSkylinePretty(samplingdates, sampling_hpd, type='step',
                    col=pal.dark(cpurple), fill=pal.dark(cpurple,0.7),
                    xaxis=TRUE, yaxis=TRUE, axispadding=0, ylims=c(0,samplingmax),
                    ylab="Sampling proportion", xlab="", xline=2, side=2, yline=2.5,
                    lwd=1, cex.axis=0.7, cex.label=1.0, xticks=xticks, xticklabels=NA,
                    new=FALSE, add=TRUE)
  
  # Box around plot
  rect(min(xticks), 0, max(xticks), samplingmax, xpd=TRUE)
  
  # Plot X-axis
  text(x=xticks[1:(length(xticks)-1)], y=-0.02*samplingmax, monthlabs[1:(length(monthlabs)-1)], xpd=TRUE, srt=0, pos=1, cex=0.7)
  text(x=xticks[yearstarts], y=-0.1*samplingmax, yearlabs, xpd=TRUE, srt=0, pos=1, cex=0.7)
  
  # Add empirical sampling proportion 
  if (!is.na(samplingfile)) {
    empiricalsampling <- data.frame(t(read.csv(samplingfile, row.names = 1)))
    dates <- getDayDate(empiricalsampling$yeardate)
    lines(c(dates[1], dates, enddate), c(0, empiricalsampling$combined, empiricalsampling$combined[nrow(empiricalsampling)]), type='s',col=pal.dark(cred),lwd=1, lty=2)
  }
  
  dev.off()
}