# Bayesian Modeling of Interaction between Features in Sparse Multivariate Count Data with Application toMicrobiome Study

ZI-MLN (Zero Inflated Multivariate rounded Log Normal Model) is a Bayesian model for the analysis of next-generation sequencing microbiome abundance data. It models interaction between microbial features in a community directly using observed count data, in the presence of covariates and excess zeros.

For more information, read the paper: ... (link)

## Installation

### Install R

ZI-MLN requires R 3.6 or greater to reproduce the tables and graphics from the paper and to compare the performance of ZI-MLN and its comparators implemented by `edgeR` and `metagenomeSeq`.

Download and install R from [https://www.r-project.org/](https://www.r-project.org/).

Once installed, open R from the terminal with `R` and run the following command:

```
install.packages(c("statmod", "GIGrvg", "extraDistr", "mvtnorm"))
```

Rstudio is a more user-friendly platform to implement R code. Download from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/).


## Simulation study

Run `without covariate.R` to simulate the model when there are no covariates on artificially-generated data. Sample size and OTU columns can be specified at the first through n and J. Other parameters specification are as discussed in the simulation part of the paper. 

Run `with covariate.R` to simulate the model on artificially-generated data. A single simulation replicate takes around 0.6 hours (without) / 1.1 hours (with) to run on a single core of a 2.6 GHz Intel Core i7 processor. 

## Data analysis

Run `with covariate.R` to fit a model to microbiome data from a designed experiment. A call to this script needs the data directory containing two RData:

- `X.RData`: related covariates in `n` samples (rows) by `p` covariates (columns)
- `Y.RData`: microbiome OTU data in `n` samples (rows) by `J` OTUs (columns)

For comparators, ZI-MLN without $\Lambda$ is just removing the latent factor part. The simplified model takes less tie to implement and includes fewer parameters to estimate. 
