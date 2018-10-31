# Plots the skylines from the raw log files
# (before truncation or combining chains) 
#
# This is just to get an idea of what the results look like and not meant to be final!
# (no idea if runs converged or not)
#

rm(list = ls())
library(bdskytools)
library(yaml)
source("SkylineUtils.R")

dateorigin <- "1970-01-01"

datapath     <- "../data/"
samplingpath <- "../results/datasets/proportions/"
logbase      <- "../results/euler/"
configbase   <- "../results/config/"
outputbase   <- "../results/figures/preview/"

plotstart <- as.Date("2013-12-01")

models <- list("BDSKYSA.BMT.relaxedclock",
               "BDSKYSA.BMT.relaxedclock.rootcondition",
               "BDSKYSA.HKY+G+F.relaxedclock",
               "BDSKYSA.HKY+G+F.relaxedclock.rootcondition",
               "BDSKYSA.HKY+G+F.relaxedclock.treeslicer")

# This is dirty and dangerous and dates should really be read from sequence files linked from yaml config files
lastdates <- list("EBOV-ALL"   =as.Date("2015-10-24"),
                  "EBOV-SUB"   =as.Date("2015-10-24"),
                  "EBOV-SUBBIG"=as.Date("2015-10-24"),
                  "EBOV-EARLY" =as.Date("2014-08-03"),
                  "SLE-E"      =as.Date("2015-03-10"),
                  "SLE-W"      =as.Date("2015-08-06"),
                  "SLE-SE"     =as.Date("2015-03-10"),
                  "LBR"        =as.Date("2015-02-14"))

monthlysampling <- list("EBOV-ALL"   ="all-monthly.csv",
                        "EBOV-SUB"   ="special-suball-monthly.csv",
                        "EBOV-SUBBIG"="special-subbig-monthly.csv",
                        "EBOV-EARLY" ="special-early-monthly.csv",
                        "SLE-E"      ="sle-east-monthly.csv",
                        "SLE-W"      ="sle-west-monthly.csv",
                        "SLE-SE"     ="sle-southeast-monthly.csv",
                        "LBR"        ="lbr-monthly.csv")

binned10sampling <- list("EBOV-ALL"   ="all-10bins.csv",
                         "EBOV-SUB"   ="special-suball-10bins.csv",
                         "EBOV-SUBBIG"="special-subbig-10bins.csv",
                         "EBOV-EARLY" ="special-early-10bins.csv",
                         "SLE-E"      ="sle-east-10bins.csv",
                         "SLE-W"      ="sle-west-10bins.csv",
                         "SLE-SE"     ="sle-southeast-10bins.csv",
                         "LBR"        ="lbr-10bins.csv")

reported <- list("SLE-E"      ="sierraleone_eastern_cases.csv",
                 "SLE-W"      ="sierraleone_western_cases.csv",
                 "SLE-SE"     ="sierraleone_southeast_cases.csv",
                 "LBR"        ="liberia_cases.csv")


########################################################################################################

sle_cases <- read.table(paste0(datapath,"cases/sierraleone_cases.csv"), header=TRUE)
gin_cases <- read.table(paste0(datapath,"cases/guinea_cases.csv"), header=TRUE)
lbr_cases <- read.table(paste0(datapath,"cases/liberia_cases.csv"), header=TRUE)
all_cases <- sle_cases + gin_cases + lbr_cases

for (model in models) {
  
  configpath <- paste0(configbase, model)
  for (configfile in list.files(configpath, pattern="*.cfg")) {
    parts <- strsplit(configfile, "\\.")[[1]]
    dataset <- parts[1]
    reproductivenumberchanges <- parts[length(parts)-2]
    samplingproportionchanges <- parts[length(parts)-1]
    
    logpath   <- paste0(logbase,model,".",reproductivenumberchanges,".",samplingproportionchanges)
    logprefix <- paste(dataset,model,reproductivenumberchanges, samplingproportionchanges, sep=".")
    for (logfile in list.files(logpath, pattern="*.log")) {
      if (grepl(logprefix, logfile, fixed=TRUE)) {
        seed <- gsub("\\.log", "", strsplit(logfile, "_")[[1]][2])
        
        print(paste("Plotting",logfile))
        dir.create(paste0(outputbase,dataset), recursive=TRUE, showWarnings=FALSE)
        
        if (samplingproportionchanges == "SM") {
          samplingfile <- paste0(samplingpath,monthlysampling[[dataset]])
        } else
          if (samplingproportionchanges == "S10") {
            samplingfile <- paste0(samplingpath,binned10sampling[[dataset]])
          } else
            samplingfile <- NA
        #ifelse((samplingproportionchanges == "SM" || samplingproportionchanges == "S10"), paste0(samplingpath, monthlysampling[[dataset]]), NA)
        
        if (dataset %in% names(reported)) {
          cases <- read.table(paste0(datapath,"cases/", reported[[dataset]]), header=TRUE)
        } else {
          cases <- all_cases
        }
        
        plotSkyline(configfile  =paste(configpath, configfile, sep="/"), 
                    logfile     =paste(logpath, logfile, sep="/"), 
                    samplingfile=samplingfile,
                    outfile     =paste0(outputbase, dataset, "/", logprefix,"_",seed,".pdf"),
                    reported    =cases,
                    startdate   =plotstart,
                    enddate     =lastdates[[dataset]])            
        
      }
    }
  }
  
}