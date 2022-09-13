# Zero Inflated Multivariate rounded Log Normal Model

ZI-MLN (Zero Inflated Multivariate rounded Log Normal Model) is a Bayesian model for the analysis of next-generation sequencing microbiome abundance data. It models interaction between microbial features in a community directly using observed count data, in the presence of covariates and excess zeros.

For more information, read the paper: Bayesian Modeling of Interaction between Features in Sparse Multivariate Count Data with Application to Microbiome Study (to appear in AOAS)

## Installation

### Install R

ZI-MLN requires R 3.6 or greater to reproduce the tables and graphics from the paper and to compare the performance of ZI-MLN. Download and install R from [https://www.r-project.org/](https://www.r-project.org/). Once installed, open R from the terminal with R and run the following command:

```
install.packages(c("statmod", "GIGrvg", "extraDistr", "mvtnorm"))
```

Rstudio is a more user-friendly platform to implement R code. Download from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/).

## Simulation study for count table without covariates 

The implementation R code is in `without covariate.R` to analysis the data when there are no covariates on artificially-generated data. The input count table do not need normalization. And for input count table, each row is each sample and each column is each features(OTUs). Notice we have also another subject index $m=1,2...M$. For example, $m=1,2,2,3$ means that the first sample belongs to the first subject. The 2nd and 3rd sample belong to the second subject. The 4th sample belong to the third subject. A special case is considered in the paper and below is $m=1,2,3...n$, which implies each subject has their own one sample. Hyper-parameters specification are as discussed in the simulation part of the paper. 

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

After we save the simulated synthetic data and the truth parameters, we can run the ZI_MLN_without function in `without covariate.R` to analysis the data. It takes about 20 minutes on Apple Macbookpro M1Max 2021.

```
ls = ZI_MLN_without(Y, M = M, m = m)
```

To access the performance, one can use the following code to compare the posterior estimates to the truth. 

```
niter= length(ls)
burn=1:niter
library(ggplot2)
library(latex2exp) ### install if it's not in the library

############################ ri+thetaj ######################

df.riplusthetaj = matrix(NA, nrow = n*J, ncol = 4)
df.riplusthetaj.mean = df.riplusthetaj.up = df.riplusthetaj.low = df.riplusthetaj.true =  matrix(NA, nrow = n, ncol = J)

for(i in 1:n){
  for(j in 1:J){
    cccccc = sapply(burn, function(x) ls[[x]]$ri[i]+ls[[x]]$thetaj[j])
    df.riplusthetaj.low[i, j] = quantile(cccccc, probs = 0.025)
    df.riplusthetaj.mean[i, j] = mean(cccccc)
    df.riplusthetaj.up[i, j] = quantile(cccccc, probs = 0.975)
    df.riplusthetaj.true[i, j] = ri.true[i] + thetaj.true[j]
  }
}

df.riplusthetaj = c(df.riplusthetaj.low)
df.riplusthetaj = cbind(df.riplusthetaj, c(df.riplusthetaj.mean))
df.riplusthetaj = cbind(df.riplusthetaj, c(df.riplusthetaj.up))
df.riplusthetaj = cbind(df.riplusthetaj, c(df.riplusthetaj.true))
df.riplusthetaj = data.frame(df.riplusthetaj)
colnames(df.riplusthetaj) <- c('2.5%', 'Pos.Mean', '97.5%', 'True')

ggplot(df.riplusthetaj, aes(x = True)) + 
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 1) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\widehat{r_i+\\alpha_j}$"))+
  xlab(TeX("$r_i+\\alpha_j$")) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=25),
        text=element_text(size=30)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

############################## vs2/sig2 ####################

mean(sapply(burn, function(x) ls[[x]]$vs2))
mean(sapply(burn, function(x) ls[[x]]$sig2))
mean(sapply(burn, function(x) ls[[x]]$sig2 + ls[[x]]$vs2))

####################### marginal pos cor  ##################

mar.pos.cor = matrix(0, J, J)
for(ilist in burn){
  mar.pos.cor = mar.pos.cor + cov2cor(tcrossprod(ls[[ilist]]$Lambda) + diag(ls[[ilist]]$vs2+ls[[ilist]]$sig2, J))
}
mar.pos.cor = mar.pos.cor/length(burn)
true.cor = cov2cor(Omega.true+diag(vs2.true, J))
df.diff.cor = data.frame(x = mar.pos.cor[upper.tri(mar.pos.cor, diag = F)] - (true.cor)[upper.tri(true.cor, diag = F)])
ggplot(df.diff.cor, aes(x=x)) + geom_histogram(color="black", fill="white", size=1.1) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Density')+
  xlim(c(-1.5,1.5))+
  xlab(TeX("$\\hat{\\rho_{jj'}}-\\rho_{jj'}^{tr}$")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=25),
        text=element_text(size=35)) 
```

## Simulation study for count table with covariates 

First we simulate a synthentic count table with covariate. 

```
n = 70 
J = 150 
K.true = 5
zero.rate = 80/100  ### again this is sparisity in covariance matrix

seed = 3
set.seed(seed)
Lambda.true = matrix(runif(J*K.true, -3, 3), nrow = J, ncol =K.true)
for(k in 1:K.true){
  set.seed(k)
  ind = sample(1:J, size = round(J*zero.rate))
  Lambda.true[ind, k] = 0
}
sigma2.true = 1
Omega.true = Lambda.true %*% t(Lambda.true) + sigma2.true * diag(J)

ri.true = runif(n, 3, 7)
thetaj.true = runif(J, 0, 2)
vs2.true = 1
mu.true = matrix(NA, n, J)
M = 35
m = rep(1:M, 2)
true.cor <- cov2cor(Omega.true+diag(vs2.true,J))

smj.true = matrix(rnorm(M*J, 0, sqrt(vs2.true)), nrow = M, ncol = J)

X = matrix(NA, nrow =n, ncol = 2)
X[1:(n/2),1] = rnorm(n/2, 0, 1)
X[((n/2+1):n),1] = X[1:(n/2),1]
X[,2] = rep(c(1,0), each = M)
tilde.X = cbind(1, X)
X = cbind(X, rep(c(0,1), each = M))

p2 = ncol(X)
set.seed(seed)
beta.true = matrix(0, nrow = J, ncol = p2)
beta.true[,1] = rnorm(J, 0, 1)
beta.true[,2] = rnorm(J, 0, 1)
beta.true[,3] = rnorm(J, 0, 1)

for(i in 1:n){
  for(j in 1:J){
    mu.true[i, j] = ri.true[i] + thetaj.true[j] + smj.true[m[i], j] + 
      sum(X[i,] * beta.true[j,])
  }
}

y.star.true = matrix(NA, nrow = n, ncol = J)
set.seed(seed)
for(i in 1:n){
  y.star.true[i, ] = mvtnorm::rmvnorm(1, mean = mu.true[i, ], Omega.true)
}
Y = floor(exp(y.star.true))

kappa.jp.true = matrix(runif(J*3, -0.5,0), nrow =J, ncol = 3)
epsilon.ij.true = pnorm(tilde.X %*% t(kappa.jp.true))
delta.ij.true = matrix(NA, nrow = n, ncol = J)
for(i in 1:n){
  for(j in 1:J){
    delta.ij.true[i, j] = rbinom(1, 1, 1-epsilon.ij.true[i, j])
  }
}

Y = Y*delta.ij.true
```

Next we run ZI_MLN_with function in `with covariate.R` to analysis the data. It takes about 1.2hour on Apple Macbookpro M1Max 2021.

Use the following code to do posterior checking:

```
niter= length(ls)
burn=1:niter
library(ggplot2)
library(latex2exp)

######################################## vs2/sig2 ##########

mean(sapply(burn, function(x) ls[[x]]$vs2))
mean(sapply(burn, function(x) ls[[x]]$sig2))
mean(sapply(burn, function(x) ls[[x]]$sig2 + ls[[x]]$vs2))

####################### marginal pos cor  ###############################

mar.pos.cor = matrix(0, J, J)
for(ilist in burn){
  mar.pos.cor = mar.pos.cor + cov2cor(tcrossprod(ls[[ilist]]$Lambda) + diag(ls[[ilist]]$vs2+ls[[ilist]]$sig2, J))
}
mar.pos.cor = mar.pos.cor/length(burn)
true.cor = cov2cor(Omega.true+diag(vs2.true, J))
df.diff.cor = data.frame(x = mar.pos.cor[upper.tri(mar.pos.cor, diag = F)] - (true.cor)[upper.tri(true.cor, diag = F)])
ggplot(df.diff.cor, aes(x=x)) + geom_histogram(color="black", fill="white", size=1.1) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Density')+
  xlim(c(-1.5,1.5))+
  xlab(TeX("$\\hat{\\rho_{jj'}}-\\rho_{jj'}^{tr}$")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=25),
        text=element_text(size=35)) 

############################ ri ###############################

df.ri = matrix(NA, nrow = n, ncol = 4)
for(i in 1:n){
  cccccc = sapply(burn, function(x) ls[[x]]$ri[i])
  df.ri[i, 1] = quantile(cccccc, probs = 0.025)
  df.ri[i, 2] = mean(cccccc)
  df.ri[i, 3] = quantile(cccccc, probs = 0.975)
}
df.ri[,4] = ri.true
df.ri = data.frame(df.ri)
colnames(df.ri) <- c('2.5%', 'Pos.Mean', '97.5%', 'True')
ggplot(df.ri, aes(x = True)) + 
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 2) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() + labs(title = "R_i") + theme(plot.title = element_text(hjust = 0.5))

############################ thetaj ###############################

df.thetaj = matrix(NA, nrow = J, ncol = 4)

for(i in 1:J){
  cccccc = sapply(burn, function(x) ls[[x]]$thetaj[i])
  df.thetaj[i, 1] = quantile(cccccc, probs = 0.025)
  df.thetaj[i, 2] = mean(cccccc)
  df.thetaj[i, 3] = quantile(cccccc, probs = 0.975)
}

df.thetaj[,4] = thetaj.true
df.thetaj = data.frame(df.thetaj)
colnames(df.thetaj) <- c('2.5%', 'Pos.Mean', '97.5%', 'True')

ggplot(df.thetaj, aes(x = True)) + 
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 2) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() +
  # geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.01, size =0.01,
  #               position=position_dodge(0.05), colour = "red") +
  labs(title = "Theta_j") +
  theme(plot.title = element_text(hjust = 0.5))

############################ ri+thetaj ###################################

new.ls = list()
for(i in burn){
  new.ls[[i]] = matrix(ls[[i]]$ri, nrow= n, ncol = J) + matrix(ls[[i]]$thetaj, byrow = T, nrow= n, ncol = J)
}

df.rithetaj.mean = Reduce("+", new.ls) / length(new.ls)
df.rithetaj.low = apply(simplify2array(new.ls), 1:2, function(x) quantile(x,0.025))
df.rithetaj.up = apply(simplify2array(new.ls), 1:2, function(x) quantile(x,0.975))
riplusthetajtrue = matrix(ri.true, nrow= n, ncol = J) + matrix(thetaj.true, byrow = T, nrow= n, ncol = J)
df.riplusthetaj = c()
df.riplusthetaj = c(df.rithetaj.low)
df.riplusthetaj = cbind(df.riplusthetaj, c(df.rithetaj.mean))
df.riplusthetaj = cbind(df.riplusthetaj, c(df.rithetaj.up))
df.riplusthetaj = cbind(df.riplusthetaj, c(riplusthetajtrue))
df.riplusthetaj = data.frame(df.riplusthetaj)
colnames(df.riplusthetaj) <- c('2.5%', 'Pos.Mean', '97.5%', 'True')

ggplot(df.riplusthetaj, aes(x = True)) +
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 1) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\hat{r_i}+\\hat{\\theta_j}$"))+
  xlab(TeX("$r_i+\\theta_j$")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

####################### beta_jp ################################

df.beta.j3 = sapply(burn, function(x) ls[[x]]$beta[,1])
dim(df.beta.j3)
df.beta.j3.plot = apply(df.beta.j3, 1, function(x) quantile(x, 0.025))
df.beta.j3.plot = cbind(df.beta.j3.plot, rowMeans(df.beta.j3))
df.beta.j3.plot = cbind(df.beta.j3.plot, apply(df.beta.j3, 1, function(x) quantile(x, 0.975)))
df.beta.j3.plot = cbind(df.beta.j3.plot, beta.true[,1])

colnames(df.beta.j3.plot) <- c('2.5%', 'Pos.beta', '97.5%', 'True')
df.beta.j3.plot = data.frame(df.beta.j3.plot)
colnames(df.beta.j3.plot) <- c('2.5%', 'Pos.beta', '97.5%', 'True')

ggplot(df.beta.j3.plot, aes(x = True)) + 
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.1, size =0.1,
                position=position_dodge(0.05), colour = "grey") +
  geom_point(aes(y = Pos.beta, colour = "Pos beta")) + 
  geom_line(aes(x = True, y = True, colour = "True"))+
  theme_bw() +
  theme(legend.position="none") +
  xlab(TeX("$\\beta_{j3}^{tr}$")) + 
  ylab(TeX("$\\hat{\\beta_{j3}}$")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=25),
        text=element_text(size=30)) 


df.beta.j2 = sapply(burn, function(x) ls[[x]]$beta[,2]- ls[[x]]$beta [,3])
df.beta.j2.plot = apply(df.beta.j2, 1, function(x) quantile(x, 0.025))
df.beta.j2.plot = cbind(df.beta.j2.plot, rowMeans(df.beta.j2))
df.beta.j2.plot = cbind(df.beta.j2.plot, apply(df.beta.j2, 1, function(x) quantile(x, 0.975)))
df.beta.j2.plot = cbind(df.beta.j2.plot, beta.true[,2]  - beta.true[,3])
colnames(df.beta.j2.plot) <- c('2.5%', 'Pos.beta', '97.5%', 'True')
df.beta.j2.plot = data.frame(df.beta.j2.plot)
colnames(df.beta.j2.plot) <- c('2.5%', 'Pos.beta', '97.5%', 'True')
ggplot(df.beta.j2.plot, aes(x = True)) + 
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.1, size =0.1,
                position=position_dodge(0.05), colour = "grey") +
  geom_point(aes(y = Pos.beta, colour = "Pos beta")) + 
  geom_line(aes(x = True, y = True, colour = "True"))+
  theme_bw() +
  theme(legend.position="none") +
  xlab(TeX("$\\beta_{j1}^{tr}-\\beta_{j2}^{tr}$")) + 
  ylab(TeX("$\\widehat{\\beta_{j1}-\\beta_{j2}}$")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=25),
        text=element_text(size=30)) 
```

## Data analysis

Run functions like simulation cases above. Choices of some hyper-parameters  are discussed in the paper. 

$a_\phi$: controls the sparsity in the covariance matrix. We recommend not using a very small value at the first try. $a_\phi=1/2$ is a common choice and if you find it does not have enough shrinkage, just gradually decrease it. 

$a_\tau, b_\tau$: we set them to be (1,50) as default. It gives a non-informative priora and one can change to (2, 50). 

$K$: the number of sub-dimension oof covariance matrix decomposition. We set 10 as default. 

$L^r, L^\theta$: the number of mixtures in the prior of $r_i$ and $\theta_j$. We set 8, 15 as default.
