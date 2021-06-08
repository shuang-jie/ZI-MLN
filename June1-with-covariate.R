rm(list = ls())

list.of.packages <- c("factoextra", "statmod", "GIGrvg", "extraDistr", "mvtnorm", "corrplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

n = 70 
J = 150 
K.true = 5
zero.rate =50/100 # 80/100

seed = 2

################## simulation setting ###################################
set.seed(seed)
Lambda.true = matrix(runif(J*K.true, -3, 3), nrow = J, ncol =K.true)

for(k in 1:K.true){
  set.seed(k)
  ind = sample(1:J, size = round(J*zero.rate))
  Lambda.true[ind, k] = 0
}
sigma2.true = 1
Omega.true = Lambda.true %*% t(Lambda.true) + sigma2.true * diag(J)
# corrplot::corrplot(cov2cor(Omega.true), tl.pos = 'n')

library(factoextra)
res.pca = prcomp(Omega.true)
fviz_eig(res.pca, addlabels = T)
try = svd(Omega.true)
eigen(Omega.true)$values
eigen(Omega.true)$values[1] / sum(eigen(Omega.true)$values)

set.seed(seed)
ri.true = runif(n, 3, 7)
thetaj.true = runif(J, 0, 2)
vs2.true = 1
mu.true = matrix(NA, n, J)
M = 35
m = rep(1:M, 2)
true.cor <- cov2cor(Omega.true+diag(vs2.true,J))
#  corrplot::corrplot(true.cor, tl.pos = 'n')
max(true.cor[1:50, 1:50][-which(true.cor[1:50, 1:50]==1)])
min(true.cor[1:50, 1:50][-which(true.cor[1:50, 1:50]==1)])

smj.true = matrix(rnorm(M*J, 0, sqrt(vs2.true)), nrow = M, ncol = J)

X = matrix(NA, nrow =n, ncol = 2)
X[1:(n/2),1] = rnorm(n/2, 0, 1)
X[((n/2+1):n),1] = X[1:(n/2),1]
X[,2] = rep(c(1,0), each = M)
tilde.X = cbind(1, X)
X = cbind(X, rep(c(0,1), each = M))

p2 = ncol(X)
beta.p.true = matrix(rnorm(J*p2, 0, 1), nrow = J, ncol = p2)

for(i in 1:n){
  for(j in 1:J){
    mu.true[i, j] = ri.true[i] + thetaj.true[j] + smj.true[m[i], j] + 
      sum(X[i,] * beta.p.true[j,])
  }
}

y.star.true = matrix(NA, nrow = n, ncol = J)
set.seed(seed)
for(i in 1:n){
  y.star.true[i, ] = mvtnorm::rmvnorm(1, mean = mu.true[i, ], Omega.true)
}
Y = floor(exp(y.star.true))
sum(Y==0) ############ 40% zero
sum(Y==0)/(n*J)

kappa.jp.true = matrix(runif(J*3, -0.5,0), nrow =J, ncol = 3)
epsilon.ij.true = pnorm(tilde.X %*% t(kappa.jp.true))
delta.ij.true = matrix(NA, nrow = n, ncol = J)
for(i in 1:n){
  for(j in 1:J){
    delta.ij.true[i, j] = rbinom(1, 1, 1-epsilon.ij.true[i, j])
  }
}

Y = Y*delta.ij.true
sum(Y==0) ############ 40% zero
sum(Y==0)/(n*J)






#############################

K = 10

library(statmod)
library(GIGrvg)
library(extraDistr)
library(mvtnorm)
################################## MCMC setting ################################

niter = 30000
ls = list()
psi.matrix = matrix(rexp(J*K, 1/2), nrow =J, ncol = K)
# a.phi = 1/0.1/2
a.phi = 1/2

a.tau = 1
b.tau = 1/2
a.sigma = 1
b.sigma = 1
accept.target = 0.44
set.seed(seed)
tilde.phi = rgamma(J, 1, 1)
tilde.tilde.phi = log(tilde.phi)
phi.vector = tilde.phi/sum(tilde.phi)
tau.k = rgamma(K, a.tau, rate = b.tau)
sigma2 = rgamma(1, a.sigma, b.sigma)
eta = rmvnorm(n, rep(0, K), diag(1, K))
Lambda = matrix(0, nrow = J, ncol =K)
sigma.proposal = rep(sqrt(2.38^2/J), J)
naccept = true.naccept = rep(0, J)
accept = rep(0, J)
delta = 1


ql = seq(0.5, 0.99,by = 0.01)
ql.matrix = sapply(1:n, function(x) quantile(Y[x,], ql))
dl.new = c()
for(i in 1:nrow(ql.matrix)){
  ql.medj = median(ql.matrix[i,])
  dl = median(abs(ql.matrix[i,] - ql.medj))
  dl.new[i] = dl
}
find.best.l = c()
for(i in 1:(length(dl.new)-1)){
  if(dl.new[i+1]> 1.1*dl.new[i]){
    find.best.l = cbind(find.best.l, i)
  }
}
best.l = ql[find.best.l[1]]
hat.ri = c()
for(i in 1:n){
  hat.ri[i] =log(0.01+sum(Y[i,Y[i,] < quantile(Y[i,], best.l)]))
}

hat.thetaj = (log(colSums(Y+0.01)) - sum(hat.ri))/n
ls = list()
## hyper-parameter
vr = mean(hat.ri); Lr = 8; a.phi.r = 3; a.omega = b.omega = 5; #ur2 = var(hat.ri)
ur2 = 1
a.xi = vr
sigma2.xi = 1

vr.theta = mean(hat.thetaj); Lr.theta = 10; a.phi.theta = 3; a.omega.theta = b.omega.theta = 5; #ur2.theta = var(hat.thetaj)
ur2.theta = 1
a.xi.theta = vr.theta
sigma2.xi.theta = 1

############################## starting point ###########################
ri = hat.ri
Si1 = sample(1:Lr, n, replace = T)
Si2 = sample(0:1, n, replace = T)
xi = rep(0, Lr)
omega.l = rbeta(Lr, a.omega, b.omega)

thetaj = hat.thetaj
Sj1 = sample(1:Lr.theta, J, replace = T)
Sj2 = sample(0:1, J, replace = T)
xi.theta = rep(0, Lr.theta)
omega.l.theta = rbeta(Lr.theta, a.omega.theta, b.omega.theta)
y.star = matrix(0, nrow = n, ncol = J)

epsilon.ij = matrix(runif(n*J), nrow = n, ncol =J)
delta.matrix = matrix(1, nrow = n, ncol = J)

posterior.phi = function(t.phi){
  # caluclate the posterior likelihood
  or.phi = t.phi/sum(t.phi)
  tt.phi = log(t.phi)
  log.prior = sum(dgamma(t.phi, a.phi, 1, log = T))
  log.like = sum(dlaplace(c(Lambda), mu = 0, 
                          sigma = rep(or.phi, K)*rep(tau.k, each = J), log=T))
  log.pos = log.prior + log.like
  Jo = t.phi[j]
  pos.final = log.pos + log(Jo)
  return(pos.final)
}

p1 = ncol(X)
p2 = ncol(tilde.X)-1

mu_beta = matrix(rep(0, p1), ncol = 1); Sigma_beta = diag(3, p1)

Z.matrix = matrix(0, nrow = n, ncol = J)

kappa.matrix = matrix(0, nrow = J, ncol = p2+1)
mu_kappa = matrix(rep(0, p2+1), ncol = 1); Sigma_kappa = diag(3, p2+1)

beta.p =  matrix(0, nrow = J, ncol =p1)
smj.matrix = matrix(0, nrow = M, ncol =J)
vs2 = 1
a.vs2 = 3
b.vs2 = 3

sij.matrix = matrix(0, nrow= n, ncol = J)
for(i in 1:n){
  for(j in 1:J){
    sij.matrix[i, j] = smj.matrix[m[i], j]
  }
}

######################## run code ##############################################

for(i in 1:niter){
  
  
  ### Update delta_ij
  for(isample in 1:n){
    for(j in 1:J){
      if( Y [isample, j] >0){
        delta.matrix[isample, j] = 1
      } else{
        prob.0 = epsilon.ij[isample, j]
        prob.mean = ri[isample] + thetaj[j] + sum(Lambda[j,] * eta[isample,]) + sum(X[isample,] * beta.p[j,]) + sij.matrix[isample, j]
        prob.1 = (1-prob.0) * pnorm(0, mean = prob.mean, sd = sqrt(sigma2))
        norm.prob.1 = prob.1 / (prob.0 + prob.1)
        delta.matrix[isample, j] = rbinom(1, 1, prob=norm.prob.1) 
      }
    }
  }
  
  ### Z.matrix
  for(isample in 1:n){
    for(j in 1: J){
      if(delta.matrix[isample, j] == 1){
        Z.matrix[isample, j] =  truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = tilde.X[isample, ] %*% kappa.matrix[j,], sd = 1)
      } else {
        Z.matrix[isample, j] = truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = tilde.X[isample, ] %*% kappa.matrix[j,], sd = 1)
      }
    }
  }
  
  ### update kappa 
  pos.var.kappa = solve( t(tilde.X) %*% tilde.X + Sigma_kappa)
  for(j in 1:J){
    kappa.matrix[j, ] = mvtnorm::rmvnorm(1, mean =  
                                           pos.var.kappa %*% (solve(Sigma_kappa) %*% mu_kappa + t(tilde.X) %*% Z.matrix[,j,drop=F] ), sigma = pos.var.kappa)
  }
  
  ### comput epsilon.ij
  epsilon.ij = pnorm(tilde.X %*% t(kappa.matrix))
  
  ### update smj
  for(msample in 1:M){
    for(j in 1:J){
      isample = which(m == msample)
      if(length(which(delta.matrix[isample, j] == 1)) == 0){
        smj.matrix[msample, j] = rnorm(1, 0, sqrt(vs2))
      } else{
        isample = isample[which(delta.matrix[isample, j] == 1)]
        smj.var = 1 / (1/vs2 + length(isample)/sigma2)
        smj.mean = smj.var * ( sum(y.star[isample, j]- thetaj[j]) - sum(ri[isample])  - sum(Lambda[j,] * colSums(eta[isample,, drop =F]))
                               - sum(X[isample, ] %*% beta.p[j,]) )  /sigma2
        smj.matrix[msample, j] =  rnorm(1, smj.mean, sqrt(smj.var))
      }
    }
  }
  
  sij.matrix = matrix(0, nrow = n, ncol = J)
  for(isample in 1:n){
    for(j in 1:J){
      sij.matrix[isample, j] = smj.matrix[m[isample], j]
    }
  }
  
  vs2 = 1/rgamma(1, a.vs2+M*J/2, b.vs2+sum(smj.matrix^2)/2)
  
  ### impute continuous latent variable
  for(isample in 1:n){
    for(j in 1:J){
      if(delta.matrix[isample, j]!=0){
        mean_uij = ri[isample] + thetaj[j]+ smj.matrix[m[isample], j] + sum(Lambda[j,] * eta[isample,])+ 
          X[isample, ] %*% beta.p[j,]
        uij = runif(1, min = pnorm(log(Y[isample, j]), mean = mean_uij, sd=sqrt(sigma2)),
                    max = pnorm(log(Y[isample, j]+1), mean = mean_uij, sd=sqrt(sigma2)))
        if(uij == 1) {uij = uij - 9.9*10^-16}
        if(uij == 0) {uij = 9.9*10^-16}
        y.star[isample,j] = qnorm(uij, mean = mean_uij, sd=sqrt(sigma2))
      }
    }
  }
  
  ### update Omega related
  
  ### Lambda
  for(j in 1:J){
    index.lambda = which(delta.matrix[,j]!=0)
    Dj = diag(psi.matrix[j, ] * phi.vector[j]^2 * tau.k^2)
    
    if(inherits(try(solve(solve(Dj)+1/sigma2 * t(Lambda[index.lambda,, drop=F]) %*% Lambda[index.lambda,, drop=F]), silent = T), "try-error")){
      Sigmaj = Dj - Dj %*% solve(Dj + sigma2 * solve( t(eta[index.lambda,, drop=F]) %*% eta[index.lambda,, drop=F])) %*% Dj
    } else {
      Sigmaj = solve(solve(Dj)+1/sigma2 * t(eta[index.lambda,, drop=F]) %*% eta[index.lambda,, drop=F])
    }
    
    muj = Sigmaj %*% t(eta[index.lambda,, drop=F]) %*% matrix(y.star[index.lambda,j]- ri[index.lambda] - thetaj[j] - smj.matrix[m[index.lambda],j] - 
                                                                (X %*% beta.p[j,])[index.lambda,], ncol=1 ) /sigma2
    if(isSymmetric(Sigmaj) == F){Sigmaj[upper.tri(Sigmaj)] <- t(Sigmaj)[upper.tri(Sigmaj)]}
    Lambda[j,] = rmvnorm(1, muj, Sigmaj)
  }
  
  ### Sigma2
  inter.RSS = y.star - sij.matrix - matrix(rep(ri,J), ncol = J) -  X %*% t(beta.p) - 
    matrix(rep(thetaj, each = n), ncol = J) - eta %*% t(Lambda)
  inter.RSS[which(delta.matrix == 0, arr.ind = T)] = 0
  sigma2 = 1/rgamma(1, a.sigma+sum(delta.matrix)/2, rate = b.sigma + 1/2 * sum(inter.RSS^2))
  
  ### etai
  for(isample in 1:n){
    index.etai = which(delta.matrix[isample, ]!=0)
    Dj = diag(1, nrow = K)
    if(inherits(try(solve(diag(1, nrow = K) + 1/sigma2 * t(Lambda[index.etai,]) %*% Lambda[index.etai,]), silent = T), "try-error")){
      Sigma.etai = Dj - Dj %*% solve(Dj + sigma2 * solve( t(Lambda[index.etai,]) %*% Lambda[index.etai,])) %*% Dj
    } else {
      Sigma.etai = solve(diag(1, nrow = K) + 1/sigma2 * t(Lambda[index.etai,]) %*% Lambda[index.etai,])
    }
    
    #Sigma.etai = solve(diag(1, nrow = K) + 1/sigma2 * t(Lambda[index.etai,]) %*% Lambda[index.etai,])
    mu.etai = Sigma.etai %*% t(Lambda[index.etai,]) %*% matrix(y.star[isample,index.etai] - ri[isample] - thetaj[index.etai] - sij.matrix[isample, index.etai]
                                                               -  (beta.p %*% X[isample,])[index.etai], ncol=1) /sigma2 
    if(isSymmetric(Sigma.etai) == F){Sigma.etai[upper.tri(Sigma.etai)] <- t(Sigma.etai)[upper.tri(Sigma.etai)]}
    eta[isample, ] = rmvnorm(1, mean = mu.etai, Sigma.etai)
  }
  
  ### psi
  for(j in 1:J){
    for(k in 1:K){
      psi.matrix[j,k] = 1/ rinvgauss(1, mean = phi.vector[j]*tau.k[k] / abs(Lambda[j,k]), 
                                     dispersion=1)
    }
  }
  
  ### tau.k
  for(k in 1:k){
    chi = 2* sum(abs(Lambda[,k]) / phi.vector)
    tau.k[k] = rgig(1, lambda = a.tau-J, psi = 2*b.tau, chi = chi)
  }
  
  ### M-H: phi.j
  delta = min(0.01, 1/sqrt(i))
  for(j in 1:J){
    
    tilde.tilde.phi.proposal = tilde.tilde.phi
    tilde.tilde.phi.proposal[j] = rnorm(1, mean = tilde.tilde.phi[j], sd = sigma.proposal[j]) 
    tilde.phi.proposal = exp(tilde.tilde.phi.proposal)
    phi.proposal = tilde.phi.proposal/sum(tilde.phi.proposal)
    diff.log.pos = posterior.phi(tilde.phi.proposal) - posterior.phi(tilde.phi)
    accept[j] <- min(1, exp(diff.log.pos))
    
    if(runif(1) < accept[j]){
      tilde.tilde.phi = tilde.tilde.phi.proposal
      tilde.phi = tilde.phi.proposal
      phi.vector = phi.proposal
      naccept[j] = naccept[j] +1
      true.naccept[j] = true.naccept[j] +1
    }
    if(i %% 50 == 0){
      ind <- ifelse(naccept[j]/50>accept.target, 1, -1)
      sigma.proposal[j] = sigma.proposal[j] * exp(delta * ind)
      naccept[j] = 0
    }
  }
  
  
  ### update ri related
  for(isample in 1:n){
    index.ri = which(delta.matrix[isample, ]!=0)
    ri.prior.mean = xi[Si1[isample]] * (Si2[isample]==1) + ((vr - omega.l[Si1[isample]] * xi[Si1[isample]]) / (1-omega.l[Si1[isample]])) * (Si2[isample]==0)
    sample.mean = sum(y.star[isample,index.ri] -sij.matrix[isample,index.ri]) - sum(thetaj[index.ri]) - sum(colSums(Lambda[index.ri,, drop=F]) * eta[isample,]) -
      sum(X[isample,] * colSums(beta.p[index.ri,, drop=F]))
    ri[isample] = rnorm(1, (ri.prior.mean/ur2 + sample.mean/sigma2)/(1/ur2 + length(index.ri)/sigma2), 1/sqrt((1/ur2 + length(index.ri)/sigma2)))
  }
  
  phi.r = MCMCpack::rdirichlet(1, sapply(1:Lr, function(x) sum(Si1==x)) + a.phi.r)
  for(l in 1:Lr){
    omega.l[l] = rbeta(1, sum(Si1==l & Si2 ==1) + a.omega, 
                       sum(Si1==l & Si2 ==0) + b.omega)
    
    xi.pos.var = 1 / (1/ sigma2.xi +  sum(Si1==l & Si2 ==1)/ur2 + 
                        omega.l[l]^2 * sum(Si1==l & Si2 ==0)/(ur2 * (1-omega.l[l])^2) )
    xi.pos.mean = xi.pos.var * (a.xi/ sigma2.xi +  sum(ri[Si1==l & Si2 ==1])/ur2 +
                                  omega.l[l] * (sum(Si1==l & Si2 ==0)*vr - (1-omega.l[l]) * sum(ri[Si1==l & Si2 ==0]))/(ur2 * (1-omega.l[l])^2) )
    
    xi[l] = rnorm(1, xi.pos.mean, sqrt(xi.pos.var))
  }
  
  for(isample in 1:n){
    pil1 = dnorm(ri[isample], xi, sqrt(ur2)) * omega.l * phi.r
    pil0 = dnorm(ri[isample], (vr-omega.l*xi)/(1-omega.l), sqrt(ur2)) * (1-omega.l) * phi.r
    if(all(pil1 == 0) & all(pil0 == 0)){pil1 = pil0 = rep(1/2, Lr)}
    pi.normalized = c(pil1, pil0)/(sum(pil1+pil0))
    ind = sample(1:(2*Lr), 1, prob = pi.normalized)
    ind1 = ifelse(ind>Lr, ind-Lr, ind)
    #ind2 = ind1 %% 2
    ind2 = ifelse(ind>Lr, 0, 1)
    Si1[isample] = ind1
    Si2[isample] = ind2
  }
  
  ### update thetaj related
  for(jsample in 1:J){
    index.thetaj = which(delta.matrix[,jsample]!=0)
    thetaj.prior.mean = xi.theta[Sj1[jsample]] * (Sj2[jsample]==1) + ((vr.theta - omega.l.theta[Sj1[jsample]] * xi.theta[Sj1[jsample]]) / (1-omega.l.theta[Sj1[jsample]])) * (Sj2[jsample]==0)
    theta.sample.mean = sum(y.star[index.thetaj,jsample] - sij.matrix[index.thetaj,jsample]) - sum(ri[index.thetaj]) - 
      sum(X[index.thetaj,, drop=F] %*% beta.p[jsample,]) - sum(Lambda[jsample,] * colSums(eta[index.thetaj,,drop=F]))
    thetaj[jsample] = rnorm(1, (thetaj.prior.mean/ur2.theta + theta.sample.mean/sigma2)/(1/ur2.theta + length(index.thetaj)/sigma2), 
                            1/sqrt((1/ur2.theta + length(index.thetaj)/sigma2)))
  }
  phi.r.theta = MCMCpack::rdirichlet(1, sapply(1:Lr.theta, function(x) sum(Sj1==x)) + a.phi.theta)
  for(l in 1:Lr.theta){
    omega.l.theta[l] = rbeta(1, sum(Sj1==l & Sj2 ==1) + a.omega.theta, 
                             sum(Sj1==l & Sj2 ==0) + b.omega.theta)
    
    xi.pos.var = 1 / (1/ sigma2.xi.theta +  sum(Sj1==l & Sj2 ==1)/ur2.theta + 
                        omega.l.theta[l]^2 * sum(Sj1==l & Sj2 ==0)/(ur2.theta * (1-omega.l.theta[l])^2) )
    xi.pos.mean = xi.pos.var * (a.xi.theta/ sigma2.xi.theta +  sum(thetaj[Sj1==l & Sj2 ==1])/ur2.theta +
                                  omega.l.theta[l] * (sum(Sj1==l & Sj2 ==0)*vr.theta - (1-omega.l.theta[l]) * sum(thetaj[Sj1==l & Sj2 ==0]))/(ur2.theta * (1-omega.l.theta[l])^2) )
    
    xi.theta[l] = rnorm(1, xi.pos.mean, sqrt(xi.pos.var))
  }
  for(jsample in 1:J){
    pil1 = dnorm(thetaj[jsample], xi.theta, sqrt(ur2.theta)) * omega.l.theta * phi.r.theta
    pil0 = dnorm(thetaj[jsample], (vr.theta-omega.l.theta *xi.theta)/(1-omega.l.theta), sqrt(ur2.theta)) * (1-omega.l.theta) * phi.r.theta
    if(all(pil1 == 0) & all(pil0 == 0)){pil1 = pil0 =  rep(1/2, Lr.theta)}
    pi.normalized = c(pil1, pil0)/(sum(pil1+pil0))
    ind = sample(1:(2*Lr.theta), 1, prob = pi.normalized)
    ind1 = ifelse(ind>Lr.theta, ind-Lr.theta, ind)
    ind2 = ifelse(ind>Lr.theta, 0, 1)
    Sj1[jsample] = ind1
    Sj2[jsample] = ind2
  }
  
  ### betaj
  for(j in 1:J){
    index.betaj = which(delta.matrix[,j]!=0)
    betaj.pos.var = solve(Sigma_beta + 1/sigma2 * t(X[index.betaj,,drop=F]) %*% (X[index.betaj,,drop=F])  )
    betaj.pos.mean = betaj.pos.var %*% (solve(Sigma_beta) %*% mu_beta + 1/sigma2 * t(X[index.betaj,,drop=F]) %*%
                                          (y.star[index.betaj,j,drop=F] - ri[index.betaj] - thetaj[j] - 
                                             sij.matrix[index.betaj,j,drop=F] - eta[index.betaj, ] %*% Lambda[j,]) ) 
    beta.p[j, ] = rmvnorm(1, betaj.pos.mean, betaj.pos.var)
  }
  
  
  Omega = Lambda %*% t(Lambda) + sigma2 * diag(1, nrow = J)
  
  ls[[i]] = list(sigma2 = sigma2,Omega = Omega, Lambda = Lambda, tau.k = tau.k, 
                 phi.vector = phi.vector, eta = eta, ri = ri, thetaj = thetaj,
                 sij.matrix = sij.matrix, delta.matrix = delta.matrix, epsilon.ij = epsilon.ij, 
                 smj.matrix = smj.matrix, kappa.matrix = kappa.matrix, 
                 beta.p = beta.p, vs2 = vs2)
  
  print(i)
}
save.image("./June1-with-50-Lr8.RData")
