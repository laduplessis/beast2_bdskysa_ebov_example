

dateorigin <- "1970-01-01"

#' Plot weekly or monthly histogram of sampling dates
#' 
#' @param samplingdates   Vector of dates (need to be converted to POSIX dates)
#' @param daterange       Range of dates to plot across (optional)
#' @param ymax            Maximum of y-axis (optional)
#' @param start.on.monday If true weeks start on Mondays, otherwise on Sundays (default FALSE)
plotSamplingHist <- function(samplingdates, daterange=NULL, breaks="weeks", col=pal.dark(cblue), ymax=NA, start.on.monday=FALSE, 
                             plotcounts=TRUE, cex.months=0.8, cex.years=1.2, cex.axis=0.8, cex.counts=0.6, ...) {
  
  # Get daterange if not supplied
  if (is.null(daterange)) {
    daterange <- range(samplingdates)
  } 
  
  # Get weekly and monthly breaks across daterange
  weeks  <- as.Date(hist(daterange, breaks='weeks',  plot=FALSE, plot.on.monday=start.on.monday)$breaks, origin=dateorigin)
  months <- as.Date(hist(daterange, breaks='months', plot=FALSE)$breaks, origin=dateorigin)
  
  # Get counts in bins
  z <- hist(samplingdates, breaks=breaks, freq=TRUE, plot = FALSE, start.on.monday=start.on.monday)
  
  # Get Y limits
  if (is.na(ymax)) {
      ymax <- max(z$counts)
  }
  ylims <- range(pretty(c(0,ymax)))
  
  # Start the plot
  plot(1,type='n',xlim=range(months),ylim=ylims, xaxs='i', yaxs='i',bty='o', xaxt='n', yaxt='n', ...)
  
  # Shade months
  for (i in seq(1,length(months),by=2)) {
    rect(months[i],ylims[1], months[i+1], ylims[2], col = pal.light(cgray), border=NA)
  }
  
  # Plot axes
  axis(2, las=2, cex.axis=cex.axis)
  axis(1, at=months, labels=NA)
  
  text(x=months[1:(length(months)-1)], y=-0.04*ylims[2], format.Date(months[2:length(months)],format="%b"), xpd=TRUE, srt=0, pos=4, cex=cex.months)
  
  yearmids <- which(format.Date(months,format="%m") == "06")
  text(x=months[yearmids], y=-0.06*ylims[2], format.Date(months[yearmids],"%Y"), xpd=TRUE, srt=0, pos=1, cex=cex.years)
  
  #yearends <- which(format.Date(months,format="%m") == "12")
  #abline(v=months[yearends], lty=2)
  #axis(1,at=months[yearends],labels=NA,line=NA,tck=-0.04)
  
  # Plot histogram bins
  lines(z, border=pal.dark(cwhite), col=col)
  
  # Add counts
  if (plotcounts) {
      countlabels <- z$counts
      countlabels[which(countlabels==0)] <- NA
  
      text(z$mids, z$counts+(0.02*ylims[2]), labels=countlabels, col="red", cex=cex.counts, srt=90, adj=c(0,0.5),xpd=TRUE)
  }
  
  # Draw new box
  rect(min(months), ylims[1], max(months), ylims[2], xpd=TRUE)
  
  return(z)  
}


#' Plot proportions of metadata in a dataset
plotProportions <- function(counts, labels, col=pal.dark(cblue), title="", ymax=NA, prop=TRUE, plotcounts=TRUE, 
                            cex.axis=par("cex.axis"), cex.names=par("cex.axis"), cex.main=par("cex.main"), cex.counts=0.6, ...) { 
  
  props <- counts[labels]
  if (prop) {
    props <- props/sum(counts)
  }
  
  if (is.na(ymax))
    ymax  <- max(props, na.rm=TRUE)*1.25
  
  barplot2(props, beside=TRUE, col=col, names.arg=NULL, border="white", las=1, ylim=c(0,ymax), cex.axis=cex.axis, ...)
  if (plotcounts) {
      counts[which(counts==0)] <- NA
      #text(1:length(labels)*1.2-0.4, props[labels]+(0.02*ymax), labels=counts[labels], col="red", cex=cex.counts, srt=90, pos=3)
      text(1:length(labels)*1.2-0.5, props[labels]+(0.02*ymax), labels=counts[labels], col="red", cex=cex.counts, srt=90, adj=c(0,0.5), xpd=TRUE)
      par(xpd=FALSE)
  }
  text(1:length(labels)*1.2-0.5, -0.01*ymax, labels=labels, xpd=TRUE, cex=cex.names, adj=c(1, 0.5), srt=90)
  
  legend("topleft",  inset=c(0,0), horiz=TRUE, bty='n', xpd=TRUE, fill=NA, border=NA, title=paste0("    ",title), cex=cex.main, legend="")
}


###########################################################################################################################################
# Maps


#' Get districts and counts from a dataset (sequences or cases)
#' 
#' Also formats district names to a standard that can be printed
getDistrictCounts <- function(dataset) { 
  
  districts    <- table(dataset$location)
  
  n <- names(districts)
  for (d in names(districts)) {
    if (d == "PortLoko")
      n[which(n == d)] <- "Port Loko"
    else if (d == "WesternRural") 
      n[which(n == d)] <- "Western Rural"
    else if (d == "WesternUrban")
      n[which(n == d)] <- "Western Urban"
    else if (d == "GrandBassa")
      n[which(n == d)] <- "Grand Bassa"
    else if (d == "GrandKru")
      n[which(n == d)] <- "Grand Kru"
    else if (d == "GrandCapeMount")
      n[which(n == d)] <- "Grand Cape Mount"
    else if (d == "RiverGee")
      n[which(n == d)] <- "River Gee"
    else if (d == "RiverCess")
      n[which(n == d)] <- "River Cess"
    #else if (d == "Foya")
    #  n[which(n == d)] <- "Lofa"
    #else if (d == "Voinjama")
    #  n[which(n == d)] <- "Lofa"
  }
  names(districts) <- n
  
  ######################################################################
  # Merge Western Area equally into Western Urban and Western Rural
  if ("Western Urban" %in% names(districts) && "WesternArea" %in% names(districts))
    districts["Western Urban"] <- districts["Western Urban"] + districts["WesternArea"]/2
  if ("Western Rural" %in% names(districts) && "WesternArea" %in% names(districts))
    districts["Western Rural"] <- districts["Western Rural"] + districts["WesternArea"]/2
  
  westernarea  <- which(names(districts) == "WesternArea")
  if (length(westernarea) > 0) 
    districts <- districts[-westernarea]
  
  ######################################################################
  
  return(districts)
}

#' Format districts from map to a standard form
formatMapDistricts <- function(districts, simplify=FALSE) {
  
  for (i in 1:length(districts)) {
    districts[i] <- gsub("é","e",districts[i])
    if (districts[i] == "Conarky") districts[i] <- "Conakry"
    else if (districts[i] == "GrandBassa") districts[i] <- "Grand Bassa"
    else if (districts[i] == "GrandKru") districts[i] <- "Grand Kru"
    else if (districts[i] == "Gbapolu") districts[i] <- "Gbarpolu"
    
    if (simplify) districts[i] <- gsub(" ", "", paste(substring(districts[i],1,1), tolower(substring(districts[i],2))))
  }
  
  return(districts)
}


#' Get all level 2 districts in each level 1 region in a list data structure
#' 
#' Simplify processes to the format used by data extracted from Genbank files
getMapRegions <- function(map, simplify=TRUE) {
  
  regions <- list()
  for (region in unique(map$NAME_1)) {
    mask <- which(map$NAME_1 == region)
    districts <- map$NAME_2[mask]
    
    if (simplify) districts <- formatMapDistricts(districts, simplify)
    
    regions[region][[1]] <- districts
  }
  
  return(regions)
}


plotColorKey <- function(cols, height=1, labels=NULL, ...) {
  
  n <- length(cols)
  N <- round(n/height)
  
  cols <- c(cols, rep(NA,N-n))
  image(z=matrix(1:N, nrow=1), col=cols, xaxt="n", yaxt="n", axes=FALSE)
  
  if (!is.null(labels)) {
    y.tics <- seq(0,n,length.out=length(labels))/N
    axis(side=4, las=1, at=y.tics, labels=labels, tick=TRUE, ...)   
  }
}


plotMapIncidence <- function(districts, country, col=pal.dark(cblue), border="#FFFFFF", maxN=-1, ...) { 
  
  #cols <- colorpanel(max(maxN, max(districts)+1), low=pal.light(cgray), high=col)
  cols <- colorpanel(max(maxN, max(districts)+1), low="#E0DED9FF", high=col)
  
  mapcolors    <- rep(pal.light(cgray), length(country$NAME_2))  
  mapdistricts <- formatMapDistricts(country$NAME_2)
  
  for (district in names(districts)) {
    idx <- which(mapdistricts == district)
    #print(paste(district, idx, mapdistricts[idx],districts[district]))
    
    if (length(idx) <= 0) {
      warning(paste(district,"not found in map!"))
    } else
      mapcolors[idx] <- cols[districts[district]+1]
  }
  
  plot(country, col=mapcolors, border=border, ...)
  
  return(cols)
}


plotMapGradient <- function(counts, map, numlabels=2, roundTo10=TRUE, maxN=NA, ...) {
  
  if (is.na(maxN)) {
    n <- max(counts)  
  } else {
    n <- maxN
  }
  
  if (roundTo10) {
    n <- ceiling(n/10)*10
  } 
  
  layout(matrix(c(1,2,3), ncol=3, byrow=TRUE), widths = c(1,5,1))
  
  par(mar=c(0,0,0,0))
  plot(1,type='n',axes=FALSE, xlab="", ylab="")
  
  par(mar=c(3,0,1,0)+0.1)
  cols <- plotMapIncidence(counts, map, maxN=n, ...)
  
  keylabels=c(as.character(round(seq(0,1,length.out=numlabels)*n)))
  
  #par(mar=c(5,0,5,8)+0.1)
  par(mar=c(3,0,1,3)+0.1)
  plotColorKey(cols, height=0.25, labels=keylabels, cex.axis=0.8)
  
  layout(matrix(c(1), ncol=1, byrow=TRUE))
  
}



