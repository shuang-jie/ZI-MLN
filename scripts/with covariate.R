library(statmod)
library(GIGrvg)
library(extraDistr)
library(mvtnorm)
library(mvnfast)

### input Y: count table Y without any covariate. n by J matrix. n samples and J OTUs. 
### input X: covariate matrix. n by p matrix. usually p<n for identitity. Each row is sample and each column is covariate. 

### input sample group
M = nrow(Y)/2  ### how many different patients/ sample groups
m = rep(1:M, 2) ### sample group specification. e.g: (1,1,1, 2,2,2) first half comes from subject 1 and second half come from subject 2. 

### function

ZI_MLN_with <- function(Y, X = X, M = M, K=10, niter = 30000, seed = 3, a.phi = 1/2, a.tau = 1, b.tau = 1/50, a.sigma = 3, b.sigma = 3,
                        a.vs2 = 3, b.vs2 = 3, accept.target = 0.44, 
                        Lr = 8, a.phi.r = 3, a.omega = 5, b.omega = 5, ur2 = 1, 
                        Lr.theta = 15, a.phi.theta = 1, a.omega.theta = 5, b.omega.theta = 5, ur2.theta = 1){
  
  ls = list()
  n = nrow(Y)
  J = ncol(Y)
  tilde.X = cbind(1, X)
  psi.matrix = matrix(rexp(J*K, 1/2), nrow =J, ncol = K)
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
    hat.ri[i] =log(0.01+sum(Y[i,Y[i,] < quantile(Y[i,], best.l)] ) )
  }
  
  hat.ri = c()
  for(i in 1:n){
    hat.ri[i] =log(0.01+sum(Y[i, which(Y[i,] < quantile(Y[i,], best.l)) ] ) )
  }
  
  hat.thetaj = (log(colSums(Y+0.01)) - sum(hat.ri))/n
  vr = mean(hat.ri); a.xi = vr; sigma2.xi = 1
  vr.theta = mean(hat.thetaj); a.xi.theta = vr.theta; sigma2.xi.theta = 1
  
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
  y.star = matrix(log(Y+0.01), nrow = n, ncol = J)
  
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
  
  Z.matrix = matrix(0, nrow = n, ncol = J)
  
  kappa.matrix = matrix(0, nrow = J, ncol = p2+1)
  mu_kappa = matrix(rep(0, p2+1), ncol = 1); Sigma_kappa = diag(3, p2+1)
  
  beta.p =  matrix(0, nrow = J, ncol =p1)
  mu_beta = matrix(rep(0, p1), ncol = 1); Sigma_beta = diag(25, p1)
  smj.matrix = matrix(0, nrow = M, ncol =J)
  
  vs2 = rgamma(1, a.vs2, b.vs2)
  sij.matrix = matrix(0, nrow= n, ncol = J)
  for(i in 1:n){
    for(j in 1:J){
      sij.matrix[i, j] = smj.matrix[m[i], j]
    }
  }
  
  index.Y.0 = which(Y==0, arr.ind = T)
  index.Y.non0 = which(Y>0, arr.ind = T)
  delta.matrix[index.Y.non0] = 1
  
  for(i in 1:niter){
    
    ### Update delta_ij
    
    delta.matrix[index.Y.0] = sapply(1:nrow(index.Y.0), function(x) {
      prob.0 = epsilon.ij[index.Y.0[x, 1], index.Y.0[x, 2]]
      prob.mean = ri[index.Y.0[x, 1]] + thetaj[index.Y.0[x, 2]] + sum(Lambda[index.Y.0[x, 2],] * eta[index.Y.0[x, 1],]) + sum(X[index.Y.0[x, 1],] * beta.p[index.Y.0[x, 2],]) + sij.matrix[index.Y.0[x, 1], index.Y.0[x, 2]]
      prob.1 = (1-prob.0) * pnorm(0, mean = prob.mean, sd = sqrt(sigma2))
      norm.prob.1 = prob.1 / (prob.0 + prob.1)
      return(rbinom(1, 1, prob=norm.prob.1))
    })
    index.delta.1 = which(delta.matrix==1, arr.ind = T)
    index.delta.0 = which(delta.matrix==0, arr.ind = T)
    
    ### Z.matrix
    Z.matrix[index.delta.1] = sapply(1:nrow(index.delta.1), function(x) {
      truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = tilde.X[index.delta.1[x,1], ] %*% kappa.matrix[index.delta.1[x,2],], sd = 1)
    })
    Z.matrix[index.delta.0] = sapply(1:nrow(index.delta.0), function(x) {
      truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = tilde.X[index.delta.1[x,1], ] %*% kappa.matrix[index.delta.1[x,2],], sd = 1)
    })
    
    
    ### update kappa 
    pos.var.kappa = solve( t(tilde.X) %*% tilde.X + solve(Sigma_kappa))
    kappa.matrix = t(sapply(1:J, function(x) rmvn(1, pos.var.kappa %*% (solve(Sigma_kappa) %*% mu_kappa + t(tilde.X) %*% Z.matrix[,x,drop=F] ), pos.var.kappa)))
    
    
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
          smj.mean = smj.var * ( sum(y.star[isample, j]- thetaj[j]) - sum(ri[isample])  - sum(Lambda[j,] * colSums(eta[isample,, drop =F])) - sum(X[isample, ] %*% beta.p[j,]) )  /sigma2
          smj.matrix[msample, j] =  rnorm(1, smj.mean, sqrt(smj.var))
        }
      }
    }
    
    ### update sij.matrix
    sij.matrix = smj.matrix[m, ]
    
    ###
    vs2 = 1/rgamma(1, a.vs2+M*J/2, b.vs2+sum(smj.matrix^2)/2)
    
    ### impute continuous latent variable
    
    mean_uij = ri[index.delta.1[,1]] + thetaj[index.delta.1[,2]] + sij.matrix[index.delta.1] + (eta %*% t(Lambda))[index.delta.1] + (X %*% t(beta.p))[index.delta.1]
    u_min = pnorm(log(Y[index.delta.1]), mean = mean_uij, sd=sqrt(sigma2))
    u_max = pnorm(log(Y[index.delta.1]+1), mean = mean_uij, sd=sqrt(sigma2))
    sample_u = runif(nrow(index.delta.1), u_min, u_max)
    sample_u[sample_u==1] = 1-9.9*10^-16
    sample_u[sample_u==0] = 9.9*10^-16
    y.star[index.delta.1] =  qnorm(sample_u, mean = mean_uij, sd=sqrt(sigma2))
    
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
      
      muj = Sigmaj %*% t(eta[index.lambda,, drop=F]) %*% matrix(y.star[index.lambda,j]- ri[index.lambda] - thetaj[j] - smj.matrix[m[index.lambda],j] - (X %*% beta.p[j,])[index.lambda,], ncol=1 ) /sigma2
      if(isSymmetric(Sigmaj) == F){Sigmaj[upper.tri(Sigmaj)] <- t(Sigmaj)[upper.tri(Sigmaj)]}
      Lambda[j,] = rmvn(1, muj, Sigmaj)
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
      eta[isample, ] = rmvn(1, mu.etai, Sigma.etai)
    }
    
    ### psi
    mean.psi = matrix(phi.vector, nrow = J, ncol =K) * matrix(tau.k, nrow = J, ncol =K, byrow = T) / abs(Lambda)
    psi.matrix = matrix(1/ rinvgauss(J*K, mean = c(mean.psi), dispersion=1), nrow=J, ncol =K)
    
    ### tau.k
    
    chi.vec = colSums(2* abs(Lambda) / matrix(phi.vector, nrow = J, ncol =K))
    tau.k = sapply(1:K, function(x) rgig(1, lambda = a.tau-J, psi = 2*b.tau, chi = chi.vec[x]))
    
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
    omega.l = sapply(1:Lr, function(l) rbeta(1, sum(Si1==l & Si2 ==1)+a.omega,  sum(Si1==l & Si2 ==0) + b.omega))
    
    xi.pos.var = sapply(1:Lr, function(l) 1 / (1/ sigma2.xi +  sum(Si1==l & Si2 ==1)/ur2 + 
                                                 omega.l[l]^2 * sum(Si1==l & Si2 ==0)/(ur2 * (1-omega.l[l])^2) ))
    xi.pos.mean = sapply(1:Lr, function(l) xi.pos.var[l] * (a.xi/ sigma2.xi +  sum(ri[Si1==l & Si2 ==1])/ur2 +omega.l[l] * (sum(Si1==l & Si2 ==0)*vr - (1-omega.l[l]) * sum(ri[Si1==l & Si2 ==0]))/(ur2 * (1-omega.l[l])^2) ))
    xi = rnorm(Lr, xi.pos.mean, sqrt(xi.pos.var))
    
    for(isample in 1:n){
      pil1 = dnorm(ri[isample], xi, sqrt(ur2)) * omega.l * phi.r
      pil0 = dnorm(ri[isample], (vr-omega.l*xi)/(1-omega.l), sqrt(ur2)) * (1-omega.l) * phi.r
      if(all(pil1 == 0) & all(pil0 == 0)){pil1 = pil0 = rep(1/2, Lr)}
      pi.normalized = c(pil1, pil0)/(sum(pil1+pil0))
      ind = sample(1:(2*Lr), 1, prob = pi.normalized)
      ind1 = ifelse(ind>Lr, ind-Lr, ind)
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
    omega.l.theta = sapply(1:Lr.theta, function(l) rbeta(1, sum(Sj1==l & Sj2 ==1) + a.omega.theta, 
                                                         sum(Sj1==l & Sj2 ==0) + b.omega.theta))
    
    xi.pos.var = sapply(1:Lr.theta, function(l) 1 / (1/ sigma2.xi.theta +  sum(Sj1==l & Sj2 ==1)/ur2.theta + omega.l.theta[l]^2 * sum(Sj1==l & Sj2 ==0)/(ur2.theta * (1-omega.l.theta[l])^2) ))
    
    xi.pos.mean = sapply(1:Lr.theta, function(l) xi.pos.var[l] * (a.xi.theta/ sigma2.xi.theta +  sum(thetaj[Sj1==l & Sj2 ==1])/ur2.theta + omega.l.theta[l] * (sum(Sj1==l & Sj2 ==0)*vr.theta - (1-omega.l.theta[l]) * sum(thetaj[Sj1==l & Sj2 ==0]))/(ur2.theta * (1-omega.l.theta[l])^2) ))
    
    xi.theta = rnorm(Lr.theta, xi.pos.mean, sqrt(xi.pos.var))
    
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
      betaj.pos.var = solve(solve(Sigma_beta) + 1/sigma2 * t(X[index.betaj,,drop=F]) %*% (X[index.betaj,,drop=F])  )
      betaj.pos.mean = betaj.pos.var %*% (solve(Sigma_beta) %*% mu_beta + 1/sigma2 * t(X[index.betaj,,drop=F]) %*% (y.star[index.betaj,j,drop=F] - ri[index.betaj] - thetaj[j] - 
                                                                                                                      sij.matrix[index.betaj,j,drop=F] - eta[index.betaj, ] %*% Lambda[j,]) ) 
      beta.p[j, ] = rmvn(1, betaj.pos.mean, betaj.pos.var)
    }
    
    if(i > 15000){
      ls[[i-15000]] = list(sigma2 = sigma2, Lambda = Lambda, tau.k = tau.k, 
                           phi.vector = phi.vector, eta = eta, ri = ri, thetaj = thetaj,
                           sij.matrix = sij.matrix, delta.matrix = delta.matrix, epsilon.ij = epsilon.ij, 
                           smj.matrix = smj.matrix, kappa.matrix = kappa.matrix, 
                           beta.p = beta.p, vs2 = vs2)
      print(i)  
    }
  }
  return(ls)
}
