rm(list = ls())
source("SkylineUtils.R")

datapath     <- "../data/"
samplingpath <- "../results/datasets/proportions/"
logbase      <- "../results/xml/"
configbase   <- "../results/config/"
outputbase   <- "../results/figures/"
dir.create(outputbase, recursive = TRUE, showWarnings = FALSE)

plotstart <- as.Date("2013-12-01")

lastdates <- list("EBOV-ALL"   =as.Date("2015-10-24"),
                  "EBOV-SUB"   =as.Date("2015-10-24"),
                  "EBOV-SUBBIG"=as.Date("2015-10-24"),
                  "SLE-E"      =as.Date("2015-03-10"),
                  "SLE-SE"     =as.Date("2015-03-10"),
                  "LBR"        =as.Date("2015-02-14"))

monthlysampling <- list("EBOV-ALL"   ="all-monthly.csv",
                        "EBOV-SUB"   ="special-suball-monthly.csv",
                        "EBOV-SUBBIG"="special-subbig-monthly.csv",
                        "SLE-E"      ="sle-east-monthly.csv",
                        "SLE-SE"     ="sle-southeast-monthly.csv",
                        "LBR"        ="lbr-monthly.csv")

binned10sampling <- list("EBOV-ALL"   ="all-10bins.csv",
                         "EBOV-SUB"   ="special-suball-10bins.csv",
                         "EBOV-SUBBIG"="special-subbig-10bins.csv",
                         "SLE-E"      ="sle-east-10bins.csv",
                         "SLE-SE"     ="sle-southeast-10bins.csv",
                         "LBR"        ="lbr-10bins.csv")

reported <- list("SLE-E"      ="sierraleone_eastern_cases.csv",
                 "SLE-SE"     ="sierraleone_southeast_cases.csv",
                 "LBR"        ="liberia_cases.csv")

marks <- as.Date(c("2013-12-26",        # Onset of symptoms in index case
                   #"2013-12-28",        # Death of index case
                   "2014-03-23",        # WHO declares Ebola outbreak
                   #"2014-03-30",        # First cases diagnosed in LBR
                   #"2014-05-25",        # First cases diagnosed in SLE
                   #"2014-06-11",        # SLE border closures
                   #"2014-06-21",        # MSF declares 2nd wave totally out of control
                   #"2014-07-27",        # LBR border closures
                   #"2014-07-30",        # Quarantines in LBR, military deployed in SLE
                   "2014-08-08",        # WHO declares Public Health Emergency of International Concern
                   #"2014-08-09",        # GIN border closures
                   #"2014-09-22",        # 150-bed treatment centre opened in Monrovia
                   "2015-05-09"         # LBR first declared Ebola free
                   #"2015-06-29",        # New cluster of cases in LBR
                   #"2015-11-07",        # SLE first declared Ebola free
                   #"2015-12-29"         # GIN first declared Ebola free
))



#######################################################################################################################
# R20.SM
datasets <- c("EBOV-SUBBIG")
model    <- "BDSKYSA.BMT.relaxedclock"
run      <- "R20.SM"

sle_cases <- read.table(paste0(datapath,"cases/sierraleone_cases.csv"), header=TRUE)
gin_cases <- read.table(paste0(datapath,"cases/guinea_cases.csv"), header=TRUE)
lbr_cases <- read.table(paste0(datapath,"cases/liberia_cases.csv"), header=TRUE)
all_cases <- sle_cases + gin_cases + lbr_cases

configpath <- paste0(configbase, model, "/")
for (dataset in datasets) {
  
    logpath  <- paste0(logbase, "/", model, ".", run,"/")
    filebase <- paste(dataset, model, run, sep=".")
    
    plotSkyline(configfile=paste0(configpath, filebase, ".cfg"), 
                logfile   =paste0(logpath, filebase,".combined_25_100000.log"), 
                samplingfile =paste0(samplingpath, monthlysampling[[dataset]]),
                outfile   =paste0(outputbase, filebase),
                reported  =all_cases,
                startdate =plotstart,
                enddate   =lastdates[[dataset]], 
                plotMarks = marks, plotMRCA=FALSE, plotOrigin=TRUE)            
    
}



#######################################################################################################################
# R10.SM
datasets <- c("EBOV-SUBBIG", "EBOV-SUB")
model    <- "BDSKYSA.BMT.relaxedclock"
run      <- "R10.SM"

sle_cases <- read.table(paste0(datapath,"cases/sierraleone_cases.csv"), header=TRUE)
gin_cases <- read.table(paste0(datapath,"cases/guinea_cases.csv"), header=TRUE)
lbr_cases <- read.table(paste0(datapath,"cases/liberia_cases.csv"), header=TRUE)
all_cases <- sle_cases + gin_cases + lbr_cases

configpath <- paste0(configbase, model, "/")
for (dataset in datasets) {
  
  logpath  <- paste0(logbase, "/", model, ".", run,"/")
  filebase <- paste(dataset, model, run, sep=".")
  
  plotSkyline(configfile=paste0(configpath, filebase, ".cfg"), 
              logfile   =paste0(logpath, filebase,".combined_25_100000.log"), 
              samplingfile =paste0(samplingpath, monthlysampling[[dataset]]),
              outfile   =paste0(outputbase, filebase),
              reported  =all_cases,
              startdate =plotstart,
              enddate   =lastdates[[dataset]], 
              plotMarks = marks, plotMRCA=FALSE, plotOrigin=TRUE)            
  
}


#######################################################################################################################
# R10.S10 (combined runs)
datasets <- c("EBOV-SUB", "LBR", "SLE-E", "SLE-SE")
model    <- "BDSKYSA.BMT.relaxedclock"
run      <- "R10.S10"

sle_cases <- read.table(paste0(datapath,"cases/sierraleone_cases.csv"), header=TRUE)
gin_cases <- read.table(paste0(datapath,"cases/guinea_cases.csv"), header=TRUE)
lbr_cases <- read.table(paste0(datapath,"cases/liberia_cases.csv"), header=TRUE)
all_cases <- sle_cases + gin_cases + lbr_cases

configpath <- paste0(configbase, model, "/")
for (dataset in datasets) {
  
  logpath  <- paste0(logbase, "/", model, ".", run,"/")
  filebase <- paste(dataset, model, run, sep=".")
  
  if (dataset %in% names(reported)) {
      cases <- read.table(paste0(datapath,"cases/", reported[[dataset]]), header=TRUE)
  } else {
      cases <- all_cases
  }
      
  plotSkyline(configfile=paste0(configpath, filebase, ".cfg"), 
              logfile   =paste0(logpath, filebase,".combined_25_100000.log"), 
              samplingfile =paste0(samplingpath, binned10sampling[[dataset]]),
              outfile   =paste0(outputbase, filebase),
              reported  =cases,
              startdate =plotstart,
              enddate   =lastdates[[dataset]], 
              plotMarks = marks, plotMRCA=FALSE, plotOrigin=TRUE)            
  
}