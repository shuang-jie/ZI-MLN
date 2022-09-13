library(statmod)
library(GIGrvg)
library(extraDistr)
library(mvtnorm)
library(mvnfast)

Y ### input Y: count table Y without any covariate. n by J matrix. n samples and J OTUs. 
X ### input X: covariate matrix. n by p matrix. usually p<n for identitity. Each row is sample and each column is covariate. 

### input sample group
M ### how many different patients/ sample groups
m ### sample group specification. e.g: (1,1,1, 2,2,2) first half comes from subject 1 and second half come from subject 2. 

### function

ZI_MLN_with <- function(Y, X = X, M = M, m = m, K=10, niter = 30000, seed = 3, a.phi = 1/2, a.tau = 1, b.tau = 1/50, a.sig = 3, b.sig = 3,
                        a.vs2 = 3, b.vs2 = 3, acc.target = 0.44, 
                        Lr = 8, a.phi.r = 3, a.w = 5, b.w = 5, ur2 = 1, 
                        L.theta = 15, a.phi.theta = 1, a.w.theta = 5, b.w.theta = 5, u2.theta = 1,
                        m_kappa = 0, sig_kappa = 3, m_beta = 0, sig_beta = 25){
  
  ls = list()
  n = nrow(Y)
  J = ncol(Y)
  zeta.m = matrix(rexp(J*K, 1/2), nrow =J, ncol = K)
  
  set.seed(seed)
  t.phi = rgamma(J, 1, 1)
  t.t.phi = log(t.phi)
  phi = t.phi/sum(t.phi)
  tau.k = rgamma(K, a.tau, b.tau)
  sig2 = rgamma(1, a.sig, b.sig)
  vs2 = rgamma(1, a.vs2, b.vs2)
  eta = rmvnorm(n, rep(0, K), diag(1, K))
  Lambda = matrix(0, nrow = J, ncol =K)
  sig.pro = rep(sqrt(2.38^2/J), J)
  nacc = true.nacc = rep(0, J)
  acc = rep(0, J)
  d = 1
  
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
  ls = list()
  ## hyper-parameter
  vr = mean(hat.ri); a.xi = vr; sig2.xi = 1
  v.theta = mean(hat.thetaj); a.xi.theta = v.theta; sig2.xi.theta = 1
  
  ri = hat.ri
  Si1 = sample(1:Lr, n, replace = T)
  Si2 = sample(0:1, n, replace = T)
  xi = rep(0, Lr)
  w.l = rbeta(Lr, a.w, b.w)
  
  thetaj = hat.thetaj
  Sj1 = sample(1:L.theta, J, replace = T)
  Sj2 = sample(0:1, J, replace = T)
  xi.theta = rep(0, L.theta)
  w.l.theta = rbeta(L.theta, a.w.theta, b.w.theta)
  y.star = matrix(log(Y+0.01), nrow = n, ncol = J)
  smj = matrix(0, nrow = M, ncol =J)
  sij = smj[m, ]
  
  eps.ij = matrix(runif(n*J), nrow = n, ncol =J)
  delta = matrix(1, nrow = n, ncol = J)

  p = ncol(X)
  tilde.X = cbind(1, X)
  
  Z = matrix(0, nrow = n, ncol = J)
  kappa = matrix(0, nrow = J, ncol = p+1)
  mu_kappa = matrix(rep(m_kappa, p+1), ncol = 1); Sigma_kappa = diag(sig_kappa, p+1)
  Sigma_kappa_inv = diag(1/sig_kappa, p+1)
  beta =  matrix(0, nrow = J, ncol =p)
  mu_beta = matrix(rep(m_beta, p), ncol = 1); Sigma_beta = diag(sig_beta, p)
  Sigma_beta_inv = diag(1/sig_beta, p)

  posterior.phi = function(t.phi) sum(dgamma(t.phi, a.phi, 1, log = T)) + sum(dlaplace(c(Lambda), mu = 0, sigma = rep(t.phi/sum(t.phi), K)*rep(tau.k, each = J), log=T)) + log(t.phi[j])
  
  index.Y.0 = which(Y==0, arr.ind = T)
  index.Y.non0 = which(Y>0, arr.ind = T)
  delta[index.Y.non0] = 1
  
  RSS = y.star - matrix(ri, n, J) - matrix(thetaj, n, J, byrow = T) - sij - tcrossprod(eta, Lambda) - tcrossprod(X, beta)
  logY = log(Y)
  logY.1 = log(Y+1)
  
  for(i in 1:niter){
    
    ### Update delta
    RSS = y.star - RSS
    delta[index.Y.0] = sapply(1:nrow(index.Y.0), function(x) {
      prob.0 = eps.ij[index.Y.0[x, 1], index.Y.0[x, 2]]
      prob.1 = (1-prob.0) * pnorm(0, mean = RSS[index.Y.0[x, 1], index.Y.0[x, 2]], sd = sqrt(sig2))
      norm.prob.1 = prob.1 / (prob.0 + prob.1)
      return(rbinom(1, 1, prob=norm.prob.1))
    })
    index.delta.1 = which(delta==1, arr.ind = T)
    index.delta.0 = which(delta==0, arr.ind = T)
    RSS = y.star - RSS
    
    ### Z
    Z[index.delta.1] = sapply(1:nrow(index.delta.1), function(x) {
      truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = tilde.X[index.delta.1[x,1], ] %*% kappa[index.delta.1[x,2],], sd = 1)
    })
    Z[index.delta.0] = sapply(1:nrow(index.delta.0), function(x) {
      truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = tilde.X[index.delta.1[x,1], ] %*% kappa[index.delta.1[x,2],], sd = 1)
    })
    
    
    ### update kappa 
    pos.var.kappa = solve( crossprod(tilde.X) + Sigma_kappa_inv)
    kappa = t(sapply(1:J, function(x) rmvn(1, pos.var.kappa %*% (Sigma_kappa_inv %*% mu_kappa + crossprod(tilde.X, Z[,x,drop=F]) ), pos.var.kappa)))
    
    ### comput eps.ij
    eps.ij = pnorm(tcrossprod(tilde.X, kappa))
    
    ### update smj
    RSS = RSS + sij
    for(msample in 1:M){
      for(j in 1:J){
        isample = which(m == msample)
        if(length(which(delta[isample, j] == 1)) == 0){
          smj[msample, j] = rnorm(1) * sqrt(vs2)
        } else{
          isample = isample[which(delta[isample, j] == 1)]
          smj.var = 1 / (1/vs2 + length(isample)/sig2)
          smj.mean = smj.var * sum(RSS[isample, j])  /sig2
          smj[msample, j] =  rnorm(1) * sqrt(smj.var) + smj.mean
        }
      }
    }
    
    ### update sij
    sij = smj[m, ]
    RSS = RSS - sij
    
    ### update vs2
    vs2 = 1/rgamma(1, a.vs2+M*J/2, b.vs2+sum(smj^2)/2)
    
    ### impute continuous latent variable
    RSS = y.star - RSS
    mean_uij = RSS[index.delta.1]
    u_min = pnorm(logY[index.delta.1], mean = mean_uij, sd=sqrt(sig2))
    u_max = pnorm(logY.1[index.delta.1], mean = mean_uij, sd=sqrt(sig2))
    sample_u = runif(nrow(index.delta.1), u_min, u_max)
    sample_u[sample_u==1] = 1-9.9*10^-16
    sample_u[sample_u==0] = 9.9*10^-16
    y.star[index.delta.1] =  qnorm(sample_u, mean = mean_uij, sd=sqrt(sig2))
    RSS = y.star - RSS
    
    ### update Omega related
    
    ### Lambda
    RSS = RSS + tcrossprod(eta, Lambda)
    for(j in 1:J){
      index.lambda = which(delta[,j]!=0)
      Q = diag(1/zeta.m[j, ] * 1/phi[j]^2 * 1/tau.k^2) + crossprod(eta[index.lambda,, drop=F]) /sig2
      U = chol(Q)
      a = crossprod(eta[index.lambda,, drop=F], RSS[index.lambda, j])/sig2
      Lambda[j,] = backsolve(U, backsolve(U, a, transpose = T) + rnorm(K))
    }
    
    ### etai
    for(isample in 1:n){
      index.etai = which(delta[isample, ]!=0)
      Q = diag(1, nrow = K) + crossprod(Lambda[index.etai,])/sig2
      a = crossprod(Lambda[index.etai,], RSS[isample, index.etai]) /sig2 
      U = chol(Q)
      eta[isample, ] = backsolve(U, backsolve(U, a, transpose = T) + rnorm(K))
    }
    RSS = RSS - tcrossprod(eta, Lambda)
    
    ### sig2
    inter.RSS = RSS
    inter.RSS[which(delta == 0, arr.ind = T)] = 0
    sig2 = 1/rgamma(1, a.sig+sum(delta)/2, rate = b.sig + 1/2 * sum(inter.RSS^2))
    
    ### psi
    mean.psi = matrix(phi, nrow = J, ncol =K) * matrix(tau.k, nrow = J, ncol =K, byrow = T) / abs(Lambda)
    zeta.m = matrix(1/ rinvgauss(J*K, mean = c(mean.psi), dispersion=1), nrow=J, ncol =K)
    
    ### tau.k
    chi.vec = colSums(2* abs(Lambda) / matrix(phi, nrow = J, ncol =K))
    tau.k = sapply(1:K, function(x) rgig(1, lambda = a.tau-J, psi = 2*b.tau, chi = chi.vec[x]))
    
    ### M-H: phi.j
    d = min(0.01, 1/sqrt(i))
    for(j in 1:J){
      t.phi.pro = t.phi
      new.pro = rnorm(1) *  sig.pro[j] + t.t.phi[j]
      t.phi.pro[j] = exp(new.pro)
      
      d.log = posterior.phi(t.phi.pro) - posterior.phi(t.phi)
      acc[j] <- min(1, exp(d.log))
      
      if(runif(1) < acc[j]){
        t.t.phi[j] = new.pro
        t.phi[j] = exp(new.pro)
        phi = t.phi/sum(t.phi)
        nacc[j] = nacc[j] +1
        true.nacc[j] = true.nacc[j] +1
      }
      if(i %% 50 == 0){
        ind <- ifelse(nacc[j]/50>acc.target, 1, -1)
        sig.pro[j] = sig.pro[j] * exp(d * ind)
        nacc[j] = 0
      }
    }
    
    ### update ri related
    RSS = RSS + matrix(ri, n, J)
    pos.v.ri = 1/(1/ur2 + rowSums(delta)/sig2)
    prior.m.ri = xi[Si1] * (Si2==1) + ((vr - w.l[Si1] * xi[Si1]) / (1-w.l[Si1])) * (Si2==0)
    sample.mean = sapply(1:n, function(x) sum(RSS[x, which(delta[x, ]!=0)]))
    ri = rnorm(n) * sqrt(pos.v.ri) + (prior.m.ri/ur2 + sample.mean/sig2) * pos.v.ri
    RSS = RSS - matrix(ri, n, J)
    
    phi.r = MCMCpack::rdirichlet(1, sapply(1:Lr, function(x) sum(Si1==x)) + a.phi.r)
    w.l = sapply(1:Lr, function(l) rbeta(1, sum(Si1==l & Si2 ==1)+a.w,  sum(Si1==l & Si2 ==0) + b.w))
    xi.pos.var = sapply(1:Lr, function(l) 1 / (1/ sig2.xi +  sum(Si1==l & Si2 ==1)/ur2 + w.l[l]^2 * sum(Si1==l & Si2 ==0)/(ur2 * (1-w.l[l])^2) ))
    xi.pos.mean = sapply(1:Lr, function(l) xi.pos.var[l] * (a.xi/ sig2.xi +  sum(ri[Si1==l & Si2 ==1])/ur2 +w.l[l] * (sum(Si1==l & Si2 ==0)*vr - (1-w.l[l]) * sum(ri[Si1==l & Si2 ==0]))/(ur2 * (1-w.l[l])^2) ))
    xi = rnorm(Lr) * sqrt(xi.pos.var) + xi.pos.mean
    
    for(isample in 1:n){
      pil1 = dnorm(ri[isample], xi, sqrt(ur2)) * w.l * phi.r
      pil0 = dnorm(ri[isample], (vr-w.l*xi)/(1-w.l), sqrt(ur2)) * (1-w.l) * phi.r
      pi.normalized = c(pil1, pil0)/(sum(pil1+pil0))
      ind = sample(1:(2*Lr), 1, prob = pi.normalized)
      Si2[isample] = ifelse(ind>Lr, 0, 1)
      Si1[isample] = ifelse(ind>Lr, ind-Lr, ind)
    }
    
    ### update thetaj related
    RSS = RSS + matrix(thetaj, n, J, byrow = T)
    pos.v.theta = 1/(1/u2.theta + colSums(delta)/sig2)
    thetaj.prior.mean = xi.theta[Sj1] * (Sj2==1) + ((v.theta - w.l.theta[Sj1] * xi.theta[Sj1]) / (1-w.l.theta[Sj1])) * (Sj2==0)
    theta.sample.mean = sapply(1:J, function(x) sum(RSS[which(delta[,x]!=0), x]))
    thetaj = rnorm(J) * sqrt(pos.v.theta) + (thetaj.prior.mean/u2.theta + theta.sample.mean/sig2) * pos.v.theta
    RSS = RSS - matrix(thetaj, n, J, byrow = T)
    
    phi.r.theta = MCMCpack::rdirichlet(1, sapply(1:L.theta, function(x) sum(Sj1==x)) + a.phi.theta)
    w.l.theta = sapply(1:L.theta, function(l) rbeta(1, sum(Sj1==l & Sj2 ==1) + a.w.theta, sum(Sj1==l & Sj2 ==0) + b.w.theta))
    xi.pos.var = sapply(1:L.theta, function(l) 1 / (1/ sig2.xi.theta +  sum(Sj1==l & Sj2 ==1)/u2.theta + w.l.theta[l]^2 * sum(Sj1==l & Sj2 ==0)/(u2.theta * (1-w.l.theta[l])^2) ))
    xi.pos.mean = sapply(1:L.theta, function(l) xi.pos.var[l] * (a.xi.theta/ sig2.xi.theta +  sum(thetaj[Sj1==l & Sj2 ==1])/u2.theta + w.l.theta[l] * (sum(Sj1==l & Sj2 ==0)*v.theta - (1-w.l.theta[l]) * sum(thetaj[Sj1==l & Sj2 ==0]))/(u2.theta * (1-w.l.theta[l])^2) ))
    xi.theta = rnorm(L.theta, xi.pos.mean, sqrt(xi.pos.var))
    
    for(jsample in 1:J){
      pil1 = dnorm(thetaj[jsample], xi.theta, sqrt(u2.theta)) * w.l.theta * phi.r.theta
      pil0 = dnorm(thetaj[jsample], (v.theta-w.l.theta *xi.theta)/(1-w.l.theta), sqrt(u2.theta)) * (1-w.l.theta) * phi.r.theta
      pi.normalized = c(pil1, pil0)/(sum(pil1+pil0))
      ind = sample(1:(2*L.theta), 1, prob = pi.normalized)
      Sj1[jsample] = ifelse(ind>L.theta, ind-L.theta, ind)
      Sj2[jsample] = ifelse(ind>L.theta, 0, 1)
    }
    
    ### betaj
    RSS = RSS + tcrossprod(X, beta)
    for(j in 1:J){
      index.betaj = which(delta[,j]!=0)
      Q = Sigma_beta_inv + crossprod(X[index.betaj,,drop=F]) /sig2
      U = chol(Q)
      a = Sigma_beta_inv %*% mu_beta + crossprod(X[index.betaj,,drop=F], RSS[index.betaj, j]) /sig2
      beta[j, ] = backsolve(U, backsolve(U, a, transpose = T) + rnorm(p))
    }
    RSS = RSS - tcrossprod(X, beta)
    
    if(i > floor(niter/2)){
      ls[[i-floor(niter/2)]] = list(Lambda = Lambda, sig2 = sig2, vs2 = vs2, phi = phi, tau.k = tau.k, eta = eta, ri = ri, thetaj = thetaj, sij = sij, smj = smj, beta = beta, delta = delta, eps.ij = eps.ij, kappa = kappa)
    }
  }
  return(ls)
}
