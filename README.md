# Probabilistic modelling improves relative dating from gene phylogenies

## What will you find here?
In this repository you will find the steps carried out for the Bayesian inference process designed to obtain a relative dating of evolutionary events using gene phylogenies. This process retrieves, from a set of gene trees, the distances between the specified nodes. Moreover, it normalises (approximately removes the rate) the calculated distances by dividing the raw distance by the branch lengths median of a reference clade. Finally, we posteriorly computed the posterior distribution of the Gamma parameters for the distribution of the normalised data, which we obtained that correlate to the Bayesian Molecular Clock model in [(Álvarez-Carretero, 2022)](https://doi.org/10.1038/s41586-021-04341-1).

## Structure
This process can be structured in 2 main sections, each section is explained individually inside the folder:

1. [**`01_distance_calculations`**](01_distance_calculations).
This is the module for calculating the different distances from the trees.

2. [**`02_inference`**](02_inference)
Once calculated the distances, we can perform the inference of the posterior distributions for the distances following the instructions of this section.

In this project data visualisation has been done separately. However, for obtaining the plots shown in the paper you just need the outputs from the scripts found here.

## Data
We deposited the data and the outputs of these scripts in a [Zenodo repository](www.zenodo.org).