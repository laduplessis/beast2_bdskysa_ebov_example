---
author    : Louis du Plessis
startdate : 2018/07/29
lastdate  : 2018/07/29
---

# Summary
Steps to extract sequencing data from XML file, partition into coding and noncoding regions and remove columns with too many unknowns.

## Input

1. `Makona_1610_cds_ig.xml`: Analysis file from Dudas _et al._ paper.

## Output

1. `Makona_1610_cds_ig.fas`: Sequence data from analysis file from Dudas _et al._ paper in Fasta format.
2. `Makona_1610_cds.fas`: Coding region of the alignment.
3. `Makona_1610_ig.fas`: Noncoding region of the alignment.
4. `Makona_1610_cds.trimmed.fas`: Coding region of the alignment with columns trimmed (0 removed).
5. `Makona_1610_ig.trimmed.fas`: Noncoding region of the alignment with columns with more than 95% unknowns removed (36 removed).

---

# Extract sequencing data
Raw alignment has coding and noncoding regions interspersed, without metadata about gene starts and ends.
Instead of using Genbank reference alignment use BEASTGen to extract alignment from a BEAST XML file used in Virus genomes reveal factors that spread and sustained the Ebola epidemic, Nature, 2017 (Dudas _et al._).

```
	# Run from templates directory (beastgen struggles when template not in path)

	java -jar ~/Documents/Projects/BEAST1/beast-mcmc/dist/beastgen.jar ../templates/to_fasta.template ../data/sequence/Analyses/Phylogenetic/Makona_1610_cds_ig.xml ../results/datasets/Makona_1610_cds_ig.fas

```

Split alignment manually at position **14518** (coding and uncoding) using AliView.


# Trim sites with too many unknowns

Remove all columns with more than 95% unknowns (> 1529/1610 sequences).
- None in coding region
- 36 in noncoding region

```
python msahist.py -i ../results/datasets/Makona_1610_cds.fas -o ../results/datasets/
python msahist.py -i ../results/datasets/Makona_1610_ig.fas -o ../results/datasets/

python trimsequences.py -a ../results/datasets/Makona_1610_cds.fas -H ../results/datasets/Makona_1610_cds.hist.csv -o ../results/datasets/ -c 0.95 -p Makona_1610_cds.trimmed
python trimsequences.py -a ../results/datasets/Makona_1610_ig.fas -H ../results/datasets/Makona_1610_ig.hist.csv -o ../results/datasets/ -c 0.95 -p Makona_1610_ig.trimmed

```



`msahist.py`

1. Create a histogram of characters at each site in an alignment as a `.csv` file. (Assumes a fixed alphabet of possible characters).

`trimsequences.py`

1. Remove positions in the alignment marked in a separate fasta file (usually by "X")
2. Remove all sites in the alignment with more than some cutoff of ambiguous characters (requires a histogram file).

`dropcols.py`

1. Used by `trimsequences.py` to remove sites from the alignment. Can also be used independently.
