rm(list = ls())
library(bdskytools)

logbase      <- "../results/xml/"
outputbase   <- "../results/figures/"
dir.create(outputbase, recursive = TRUE, showWarnings = FALSE)

periodlims <- c(10,20)
clocklims  <- c(0.0009, 0.0013)

#######################################################################################################################

plotParDensity <- function(posterior,  outfile, bw='sj', col=pal.dark(cblue), fill=pal.dark(cblue,0.5), ...) { 
  
  NewFig(outfile, width=1.5, aspectratio=3/2)
  par(mar=c(3,4,0,1)+0.1)
  
  par_hpd     <- getHPD(posterior)
  print(par_hpd)
  par_density <- density(posterior, bw=bw, from=par_hpd[1], to=par_hpd[3])
  
  plot(c(par_density$x), c(par_density$y), type='n', col=col, ...)
  polygon(c(par_density$x[1], par_density$x, par_density$x[length(par_density$x)]), c(0, par_density$y, 0), col=fill, border=NA)
  abline(v=par_hpd[2], col=col, lwd=1, lty=2)
  #abline(v=truth, col=pal.dark(cred), lwd=2, lty=2)
  
  dev.off()
}



#######################################################################################################################
# R20.SM
datasets <- c("EBOV-SUBBIG")
model    <- "BDSKYSA.BMT.relaxedclock"
run      <- "R20.SM"

for (dataset in datasets) {
  
  logpath  <- paste0(logbase, "/", model, ".", run,"/")
  filebase <- paste(dataset, model, run, sep=".")
  print(filebase)
  
  logfile  <- readLogfile(paste0(logpath, filebase,".combined_25_100000.log"))
  plotParDensity(365/logfile$becomeUninfectiousRate, outfile=paste0(outputbase, filebase, ".InfectedPeriod.pdf"),
                 ylab="Infected period\n(days)", xlab="", yaxt='n', col=pal.dark(cgreen), fill=pal.dark(cgreen,0.5), 
                 cex.axis=0.7, mgp=c(0.5,0.5,0), xlim=periodlims) 
  plotParDensity(logfile$rate.mean, outfile=paste0(outputbase, filebase, ".ClockRate.pdf"),
                 ylab="Clock rate (s/s/y)", xlab="", yaxt='n', col=pal.dark(cgreen), fill=pal.dark(cgreen,0.5), 
                 cex.axis=0.7, mgp=c(0.5,0.5,0), xlim=clocklims) 
  
}

#######################################################################################################################
# R10.SM
datasets <- c("EBOV-SUBBIG", "EBOV-SUB")
model    <- "BDSKYSA.BMT.relaxedclock"
run      <- "R10.SM"

for (dataset in datasets) {
  
    logpath  <- paste0(logbase, "/", model, ".", run,"/")
    filebase <- paste(dataset, model, run, sep=".")
    print(filebase)
    
    logfile  <- readLogfile(paste0(logpath, filebase,".combined_25_100000.log"))
    plotParDensity(365/logfile$becomeUninfectiousRate, outfile=paste0(outputbase, filebase, ".InfectedPeriod.pdf"),
                   ylab="Infected period\n(days)", xlab="", yaxt='n', col=pal.dark(cgreen), fill=pal.dark(cgreen,0.5), 
                   cex.axis=0.7, mgp=c(0.5,0.5,0), xlim=periodlims) 
    plotParDensity(logfile$rate.mean, outfile=paste0(outputbase, filebase, ".ClockRate.pdf"),
                   ylab="Clock rate (s/s/y)", xlab="", yaxt='n', col=pal.dark(cgreen), fill=pal.dark(cgreen,0.5), 
                   cex.axis=0.7, mgp=c(0.5,0.5,0), xlim=clocklims) 
  
}


#######################################################################################################################
# R10.S10 (combined runs)
datasets <- c("EBOV-SUB", "LBR", "SLE-E", "SLE-SE")
model    <- "BDSKYSA.BMT.relaxedclock"
run      <- "R10.S10"

for (dataset in datasets) {
  
  logpath  <- paste0(logbase, "/", model, ".", run,"/")
  filebase <- paste(dataset, model, run, sep=".")
  print(filebase)
  
  logfile  <- readLogfile(paste0(logpath, filebase,".combined_25_100000.log"))
  plotParDensity(365/logfile$becomeUninfectiousRate, outfile=paste0(outputbase, filebase, ".InfectedPeriod.pdf"),
                 ylab="Infected period\n(days)", xlab="", yaxt='n', col=pal.dark(cgreen), fill=pal.dark(cgreen,0.5), 
                 cex.axis=0.7, mgp=c(0.5,0.5,0), xlim=periodlims) 
  plotParDensity(logfile$rate.mean, outfile=paste0(outputbase, filebase, ".ClockRate.pdf"),
                 ylab="Clock rate (s/s/y)", xlab="", yaxt='n', col=pal.dark(cgreen), fill=pal.dark(cgreen,0.5), 
                 cex.axis=0.7, mgp=c(0.5,0.5,0), xlim=clocklims) 
  
}
