# BEAST v2.5 West Africa Ebola example

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1476124.svg)](https://doi.org/10.5281/zenodo.1476124)

Please contact Louis du Plessis (louis.duplessis@zoo.ox.ac.uk) for any questions.

## Summary
Phylodynamics example for BEAST 2.5 paper using data from the West African Ebola dataset. The example uses a sampled-ancestor birth-death skyline to infer population dynamics. This can also double as a bModelTest example.

Full [documentation](doc/ebov_beast2_example.pdf) is inside the `doc/` directory.

## Reproducing results
Follow the workflows in `workflows/`

1. `data_partition.md`: Extract and process coding and noncoding regions of the full alignment.
2. `extract_datasdets.rmd`: Extract particular datasets from alignment.
3. `extract_sampling.rmd`: Extract empirical sampling proportions from datasets.
4. `beast_runs.md`: Create, run and combine BEAST XMLs.
5. `plot_figures.md`: Plot output figures.


## License

The content of this project itself is licensed under the [Creative Commons Attribution 3.0 Unported license](https://creativecommons.org/licenses/by/3.0/), and the underlying source code used to format and display that content is licensed under the [GNU General Public License v3](LICENSE.md).