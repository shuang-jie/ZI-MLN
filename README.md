# Zero Inflated Multivariate rounded Log Normal Model

ZI-MLN (Zero Inflated Multivariate rounded Log Normal Model) is a Bayesian model for the analysis of next-generation sequencing microbiome abundance data. It models interaction between microbial features in a community directly using observed count data, in the presence of covariates and excess zeros.

For more information, read the paper: Bayesian Modeling of Interaction between Features in Sparse Multivariate Count Data with Application to Microbiome Study (to appear in AOAS)

## Installation

### Install R

ZI-MLN requires R 3.6 or greater to reproduce the tables and graphics from the paper and to compare the performance of ZI-MLN and its comparators implemented by `edgeR` and `metagenomeSeq`.

Download and install R from [https://www.r-project.org/](https://www.r-project.org/).

Once installed, open R from the terminal with `R` and run the following command:

```
install.packages(c("statmod", "GIGrvg", "extraDistr", "mvtnorm"))
```

Rstudio is a more user-friendly platform to implement R code. Download from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/).


## Simulation study for count table without covariates 

The implementation R code is in `without covariate.R` to analysis the data when there are no covariates on artificially-generated data. The input count table do not need normalization. And for input count table, each row is each sample and each column is each features(OTUs). Notice we have also another subject index $m$. For example, $m=1,2,2,3$ means that the first sample belongs to the first subject. The 2nd and 3rd sample belong to the second subject. The 4th sample belong to the third subject. A special case is considered in the paper and below is $m=1,2,3...n$, which implies each subject has their own one sample. Hyper-parameters specification are as discussed in the simulation part of the paper. 

A toy example is below for simulating synthentic data:


```
n = 20
J = 150
K.true = 5
zero.rate =80/100 ### notice this zero rate is for the sparsity in the covariance matrix
seed = 1
set.seed(seed)
Lambda.true = matrix(runif(J*K.true, -3, 3), nrow = J, ncol =K.true)
for(k in 1:K.true){
  set.seed(k)
  ind = sample(1:J, size = round(J*zero.rate))
  Lambda.true[ind, k] = 0
}
sig2.true = 1
Omega.true = Lambda.true %*% t(Lambda.true) + sig2.true * diag(J)
vs2.true = 1
true.cor <- cov2cor(Omega.true+diag(vs2.true,J))
ri.true = runif(n, 3, 7)
thetaj.true = runif(J, 0, 2)
mu.true = matrix(NA, n, J)
M = n
m = 1:M
sij.true = smj.true = matrix(rnorm(n*J, 0, sqrt(vs2.true)), nrow = M, ncol = J)
tilde.X = matrix(1, nrow = n, ncol = 1)
for(i in 1:n){
  for(j in 1:J){
    mu.true[i, j] = ri.true[i] + thetaj.true[j] + smj.true[m[i], j]
  }
}
y.star.true = matrix(NA, nrow = n, ncol = J)
for(i in 1:n){
  y.star.true[i, ] = mvtnorm::rmvnorm(1, mean = mu.true[i, ], Omega.true)
}
Y = floor(exp(y.star.true))
kappa.jp.true = matrix(runif(J*1, -1,0), nrow =J, ncol = 1)
eps.j.true = pnorm(tilde.X %*% t(kappa.jp.true))
delta.ij.true = matrix(NA, nrow = n, ncol = J)
for(i in 1:n){
  for(j in 1:J){
    delta.ij.true[i, j] = rbinom(1, 1, 1-eps.j.true[i, j])
  }
}
Y = Y*delta.ij.true
```


## Simulation study for count table with covariates 

Run `with covariate.R` to simulate the model on artificially-generated data. A single simulation replicate takes around 0.6 hours (without) / 1.1 hours (with) to run on a single core of a 2.6 GHz Intel Core i7 processor. 

## Data analysis

Run `with covariate.R` to fit a model to microbiome data from a designed experiment. A call to this script needs the data directory containing two RData:

- `X.RData`: related covariates in `n` samples (rows) by `p` covariates (columns)
- `Y.RData`: microbiome OTU data in `n` samples (rows) by `J` OTUs (columns)

For comparators, ZI-MLN without ![equation](https://latex.codecogs.com/gif.latex?\Lambda) is just removing the latent factor part. The simplified model takes less tie to implement and includes fewer parameters to estimate. The other two comprting methods are built in `edgeR` and `metagenomeSeq`.
