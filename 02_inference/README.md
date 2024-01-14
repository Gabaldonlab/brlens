# Data generating model inference

We designed a set of scripts which allow you to run the Bayesian
analysis that we designed for our data. First of all, the model we
assumed as the one that generates our data is the Gamma distribution.
Thus, we assume that the normalised branch lengths are distributed
following this distribution.

To infer the parameters $(\alpha, \beta)$ linked to this distribution 
we used a JAGS, a Gibbs sampler which retrieves posterior samples of the
parameters and some statistics of the distribution. We designed a series
of functions which make this process easier. We also added the Lognormal
distribution model to assess the fitting in comparison to the Gamma one.

The module is separated in 4 essential scripts:
```
03_inference/
├── exec_mcmc.R: the script to execute the analysis using custom inputs
├── mcmc_analysis.R: functions for the MCMC convergence and quality analysis
├── ML_inference.R: ML inference of the same processes
└── sampling.R: contains the Gibbs sampling functions
```

To execute the Bayesian analysis we will have to get the prepared data
which has been previously processed by running `distances_filtering.R`. 
This processing has to return a set of `RData` files which can be read 
by the file. In this `edat` is the `data.frame` containing the data for 
the events, it has the following columns:

> `seed`: the seed of the tree.
>
> `node`: in each gene tree the distance is calculated between two 
nodes: the origin node is the seed and the destination node is the one 
of interest and the one we want to date. This column refers to the 
destination node.
>
> `dist`: is the raw distance between the origin and destination nodes.
>
> `ndist`: is the normalised distance between the origin and destination nodes. For normalising the raw distance, we used the median
of branch lengths present in the primates node. This step is done when 
running [get_t2t_dist.py](../01_distance_calculations/get_t2t_dist.py) 
and [get_t2in_dist.py](../01_distance_calculations/get_t2in_dist.py).

For running the MCMC process which infers the posterior distribution
of the Gamma and Lognormal distributions you have to run the following command:

```bash
Rscript exec_mcmc.R <filtered distanes RData> <data.frame name> <grouping column> <distances column> <number of independent chains> <number of iterations> <thinning> <uniform prior maximum> <burnin proportion 0-1> <subsample proportion 0-1> <group to analyse> <output folder>
```

For instance, if we want to infer the posterior for the event in the node 44, we should run:
```
Rscript exec_mcmc.R ev_dists.RData edat node ndist 3 100000 3 100 0.1 1 44 group_44
```