---
author    : Louis du Plessis
startdate : 2018/08/29
lastdate  : 2018/10/21
---

# Summary

Steps to generate XML files for BEAST analyses and to run the files.

# Generate XML files

All models use the sampled-ancestors birth-death skyline under an uncorrelated relaxed clock model and only differ in the substitution model and the parameterisation of the birth-death skyline.


## HKY+G+F model

```
    python makeBeastXML.py -q 120 -t 4 -s 127,128,129,130 -n beastruns -i ../results/config/BDSKYSA.HKY+G+F.relaxedclock/
    python makeBeastXML.py -q 120 -t 4 -s 127,128,129,130 -n beastruns -i ../results/config/BDSKYSA.HKY+G+F.relaxedclock.rootcondition/
    python makeBeastXML.py -q 120 -t 4 -s 127,128,129,130 -n beastruns -i ../results/config/BDSKYSA.HKY+G+F.relaxedclock.treeslicer/
```

## bModelTest

```
    python makeBeastXML.py -q 120 -t 2 -s 127,128,129,130 -n beastruns -i ../results/config/BDSKYSA.BMT.relaxedclock/
    python makeBeastXML.py -q 120 -t 2 -s 127,128,129,130 -n beastruns -i ../results/config/BDSKYSA.BMT.relaxedclock.rootcondition/
```



# Run analyses on Euler

## HKY+G+F model
- Copy to Euler
- Make scripts executable `chmod +x beastruns.euler.sh`
- Run analyses: `./beastruns.euler.sh ../beast2.jar`

## bModelTest model
- Copy to euler
- Make scripts executable `chmod +x beastruns.euler.sh`
- Install bModelTest `packagemanager -add bModelTest`
- Modify `beastruns.euler.sh` to use `beast` script in `bin/` instead of the `.jar` file.
- Run analyses: `./beastruns.euler.sh ../../BEASTv2.5.0/bin/beast`

---

## Wait

_At least one week._

---

# Post-processing

## Truncate log and trees files

Truncate the two `EBOV-SUBBIG` runs that did not reach 200 million states in a week down to 100 million states (not necessary, but means all chains are the same length before combining and 100 million states are sufficient for reaching convergence with 25% burnin).

```
../../../scripts/truncate.sh EBOV-SUBBIG.BDSKYSA.BMT.relaxedclock.R20.SM 100000000
../../../scripts/truncate.sh EBOV-SUBBIG.BDSKYSA.BMT.relaxedclock.R10.SM 100000000
```


## Combine chains

When combining chains resample every 100000 and use 25% burnin on every chain.

Analyse the following runs further:

- BDSKYSA.BMT.relaxedclock.R10.S10
    - EBOV-EARLY
    - EBOV-SUB
    - LBR
    - SLE-E
    - SLE-SE
    - SLE-W
- BDSKYSA.BMT.relaxedclock.R10.SM
    - EBOV-SUB
    - EBOV-SUBBIG
- BDSKYSA.BMT.relaxedclock.R20.SM
    - EBOV-SUBBIG


## MCC trees

TreeAnnotator with 0% burn-in for combined runs (burn-in already removed), 10% burn-in for single runs (uncombined), use 0.0 posterior probability limit (default) and median heights (common ancestor heights won't work for sampled-ancestor trees).
Use low memory for bigger trees.
