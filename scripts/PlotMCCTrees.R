rm(list = ls())
library(ggtree)
library(bdskytools)


logbase      <- "../results/xml/"
outputbase   <- "../results/figures/"
dir.create(outputbase, recursive = TRUE, showWarnings = FALSE)

lastdates <- list("EBOV-ALL"   =as.Date("2015-10-24"),
                  "EBOV-SUB"   =as.Date("2015-10-24"),
                  "EBOV-SUBBIG"=as.Date("2015-10-24"),
                  "SLE-E"      =as.Date("2015-03-10"),
                  "SLE-SE"     =as.Date("2015-03-10"),
                  "LBR"        =as.Date("2015-02-14"))

#######################################################################################################################

selectFromSeq <- function(id, group, sep="\\|") {
  return(strsplit(id,split=sep)[[1]][group]) 
}

getSeqParts <- function(labels, group=2, sep="\\|") {
  return(sapply(labels, selectFromSeq, group=group, sep=sep))
}


plotTree <- function(infile, outfile, enddate, daterange) {
  
  # Read tree
  beasttree <- read.beast(infile)
  
  # Get country metadata
  seqnames <- beasttree@phylo$tip.label
  country  <- c()
  for (seq in seqnames) {
    country <- c(country, selectFromSeq(seq, 4))
  }
  country[which(country == "GIN")] <- "Guinea"
  country[which(country == "LBR")] <- "Liberia"
  country[which(country == "SLE")] <- "Sierra Leone"
  country[which(country == "NA")]  <- "Unknown"
  
  metadata <- data.frame(cbind(seqnames, country))
  colnames(metadata) <- c("Label","Country")
  
  
  # Get x-axis labels
  plotmonths <- getMonths(daterange[1], daterange[2])
  xticks     <- getYearDate(plotmonths)
  monthlabs  <- format.Date(plotmonths,format="%b\n ")
  yearstarts <- which(format.Date(plotmonths,format="%m") == "01")
  monthlabs[yearstarts] <- format.Date(plotmonths[yearstarts],format="%b\n%Y")
  
  
  # Plot tree
  p <- ggtree(beasttree,right=TRUE, mrsd=enddate, lwd=0.15)
  p <- p %<+% metadata +
    geom_tippoint(aes(color=Country), size=0.15) + 
    theme_tree2(legend.position="right", text = element_text(size=5)) + 
    scale_x_continuous(breaks=xticks, labels=monthlabs) 
  
  # Save file
  ggsave(filename=outfile, plot = p, width=3.5, height=3.5)

}



#######################################################################################################################
# R20.SM
datasets <- c("EBOV-SUBBIG")
model    <- "BDSKYSA.BMT.relaxedclock"
run      <- "R20.SM"

for (dataset in datasets) {
    
    logpath  <- paste0(logbase, "/", model, ".", run,"/")
    filebase <- paste(dataset, model, run, sep=".")
    outfile  <- paste0(outputbase, filebase,".MCC.pdf")
    
    plotTree(infile =paste0(logpath, filebase, ".combined_25_100000.MCC.tree"), 
             outfile=paste0(outputbase, filebase,".MCC.pdf"), 
             enddate=lastdates[[dataset]], 
             daterange=as.Date(c("2014-01-01","2015-12-01")))
}


#######################################################################################################################
# R10.SM
datasets <- c("EBOV-SUBBIG", "EBOV-SUB")
model    <- "BDSKYSA.BMT.relaxedclock"
run      <- "R10.SM"

for (dataset in datasets) {
  
    logpath  <- paste0(logbase, "/", model, ".", run,"/")
    filebase <- paste(dataset, model, run, sep=".")
    outfile  <- paste0(outputbase, filebase,".MCC.pdf")
    
    plotTree(infile =paste0(logpath, filebase, ".combined_25_100000.MCC.tree"), 
             outfile=paste0(outputbase, filebase,".MCC.pdf"), 
             enddate=lastdates[[dataset]], 
             daterange=as.Date(c("2014-01-01","2015-12-01")))
}


#######################################################################################################################
# R10.S10 (combined runs)
datasets <- c("EBOV-SUB", "LBR", "SLE-E", "SLE-SE")
model    <- "BDSKYSA.BMT.relaxedclock"
run      <- "R10.S10"

for (dataset in datasets) {
  
  logpath  <- paste0(logbase, "/", model, ".", run,"/")
  filebase <- paste(dataset, model, run, sep=".")
  outfile  <- paste0(outputbase, filebase,".MCC.pdf")
  
  plotTree(infile =paste0(logpath, filebase, ".combined_25_100000.MCC.tree"), 
           outfile=paste0(outputbase, filebase,".MCC.pdf"), 
           enddate=lastdates[[dataset]], 
           daterange=as.Date(c("2014-01-01","2015-12-01")))
}

