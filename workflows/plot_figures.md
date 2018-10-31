---
author    : Louis du Plessis
startdate : 2018/10/25
lastdate  : 2018/10/25
---

# Summary

Steps to plot output figures.

Scripts require `bdskytools` ([https://github.com/laduplessis/bdskytools](https://github.com/laduplessis/bdskytools)), `ggtree` ([https://bioconductor.org/packages/release/bioc/html/ggtree.html](https://bioconductor.org/packages/release/bioc/html/ggtree.html)), `yaml` (on CRAN).

# Preview figures

All runs, just to check what logfile contents look like.

## `PlotSkylinesPreview.R`

Run from `scripts/` to preview the reproductive number and sampling proportion skylines of all BEAST2 output files (for the models specified in the file). Origin time, tMRCA and empirical sampling proportion (where available) are also plotted.

The paths at the start of the file may need adjusting, depending on where output, config, empirical sampling proportion files are stored.


# Final figures

Scripts to plot figures for:

- `EBOV-SUBBIG.BDSKYSA.BMT.relaxedclock.R20.SM`
- `EBOV-SUBBIG.BDSKYSA.BMT.relaxedclock.R10.SM`
- `EBOV-SUB.BDSKYSA.BMT.relaxedclock.R10.SM`
- `EBOV-SUB.BDSKYSA.BMT.relaxedclock.R10.S10`
- `EBOV-EARLY.BDSKYSA.BMT.relaxedclock.R10.S10`
- `LBR.BDSKYSA.BMT.relaxedclock.R10.S10`
- `SLE-E.BDSKYSA.BMT.relaxedclock.R10.S10`
-	`SLE-SE.BDSKYSA.BMT.relaxedclock.R10.S10`
-	`SLE-W.BDSKYSA.BMT.relaxedclock.R10.S10`


The paths at the start of the files may need adjusting, depending on where output, config, reported cases and empirical sampling proportion files are stored.

##  Plot skyline figures: `PlotSkylines.R`

Run from `scripts/` to plot the reproductive number and sampling proportion skylines. Origin time, tMRCA, reported cases and empirical sampling proportion are also plotted.


## Plot trees: `PlotMCCTrees.R`

Run from `scripts/` to plot the MCC trees. The tip nodes are annotated by country.


## Plot posterior distributions: `PlotDensities.R`

Run from `scripts/` to plot the posterior distributions of the clock rate and the becoming uninfectious period (infected period). Distributions are cut-off at the 95% HPD limits.
