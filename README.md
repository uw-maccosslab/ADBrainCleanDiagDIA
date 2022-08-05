# ADBrainCleanDiagDIA
A Peptide-Centric Quantitative Proteomics Dataset for the Phenotypic Assessment of Alzheimer’s Disease

This package contains the datasets and code to generate the figures included in the scientific data paper:

* Level 0: Raw files
* Level 1: Skyline document grouped by batch
* Level 2: Skyline output with the integrated peak area for each peptide (row) in each replicate (column)
* Level 3A: normalized peptide abundance across all batches
* Level 3B: normalized protein abundance across all batches

## Installation
```R
# Install the latest version from GitHub
devtools::install_github("uw-maccosslab/ADBrainCleanDiagDIA")
```

## Generate figures
```R
get_fig2;
get_fig3;
get_fig4;
get_fig5a_b;
get_fig5c;
get_supfig1_2;
get_supfig3;
get_supfig5;
```
