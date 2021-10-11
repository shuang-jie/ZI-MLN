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

Run `scripts/fit-mcmc.jl` to fit a model to microbiome data from a designed experiment. A call to this script may look like the following:

```
mkdir nutnet-analysis-results
julia scripts/fit-mcmc.jl \
    --hyper nutnet-analysis/configs/hyper.yml \
    --monitor nutnet-analysis/configs/monitor-mimix.yml \
    --inits nutnet-analysis/configs/inits.yml \
    --factors 20 \
    nutnet-analysis/reduced-data \
    nutnet-analysis-results
```

The data directory must contain three files:
- `X.csv`: treatment covariates in `n` samples (rows) by `p` covariates (columns)
- `Y.csv`: microbiome abundance data in `n` samples (rows) by `K` taxa (columns)
- `Z.csv`: block identifiers in `n` samples (rows) by `q` blocking factors (columns)

Run `./nutnet-analysis/run-nutnet-analysis.sh -d nutnet-analysis/full-data -o nutnet-analysis -f 166 -i 20000 -b 10000 -t 20 -c 1` to reproduce the full NutNet data analysis, which will take aabout 26 hours on a 2.3 GHz Intel Core i5 processor.

To demo the analysis on a small dataset, replace `nutnet-analysis/full-data` with `nutnet-analysis/reduced-data` and select the number of factors (`-f`) to be 100 or fewer.

