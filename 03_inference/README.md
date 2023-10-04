# Data generating model inference

We designed a set of scripts which allow you to run the Bayesian
analysis that we designed for our data. First of all, the model we
assumed as the one that generates our data is the Gamma distribution.
Thus, we assume that the normalised branch lengths are distributed
following this distribution.

To infer the parameters $(\alpha, \beta)$ linked to this distribution 
we used a JAGS, a Gibbs sampler which retrieves posterior samples of the
parameters and some statistics of the distribution. We designed a series
of functions which make this process easier. We also added the normal
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
which has been previously processed by running `distances_filtering.R`
