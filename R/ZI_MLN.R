#' Fit the Zero-Inflated Multivariate Rounded Log-Normal (ZI-MLN) model
#'
#' Runs the Markov chain Monte Carlo (MCMC) sampler for the ZI-MLN model of
#' Zhang et al. (2023). The model infers the interaction (covariance) structure
#' between microbial features (OTUs) directly from an observed count table,
#' accommodating overdispersion, compositionality, subject-level random effects,
#' covariate effects on abundance, and excess (structural) zeros.
#'
#' A single, unified interface handles both the no-covariate and the
#' with-covariate cases: when `X = NULL` the covariate regression terms are
#' dropped and the zero-inflation component reduces to an intercept-only probit
#' model, reproducing the "without covariate" sampler.
#'
#' @param Y Integer count matrix, `n` samples (rows) by `J` OTUs (columns). Raw
#'   counts; no normalization is required. Zeros are allowed and modelled.
#' @param X Optional numeric covariate matrix, `n` by `P`, one row per sample.
#'   Used both in the mean model for abundance and (with an added intercept) in
#'   the probit model for the zero-inflation probability. `NULL` (default) fits
#'   the no-covariate model.
#' @param m Integer vector of length `n` giving the subject/group index of each
#'   sample. For example `c(1, 1, 2, 3)` means samples 1-2 come from subject 1,
#'   sample 3 from subject 2 and sample 4 from subject 3. Defaults to
#'   `seq_len(nrow(Y))` (each sample is its own subject).
#' @param M Number of distinct subjects/groups. Defaults to
#'   `length(unique(m))`.
#' @param K Factor dimension of the low-rank covariance decomposition
#'   \eqn{\Sigma = \Lambda\Lambda' + \sigma^2 I}. Fixed at a moderately large
#'   value (default 10).
#' @param niter Total number of MCMC iterations (default 20000).
#' @param burnin Number of initial iterations discarded as burn-in. Draws are
#'   stored for iterations `> burnin`. Defaults to `floor(niter / 2)`.
#' @param seed Random seed for reproducibility (default 3).
#' @param a.phi Dirichlet-Laplace concentration \eqn{a_\phi} controlling
#'   sparsity of the covariance (default 1/2). Smaller values shrink more.
#' @param a.tau,b.tau Gamma shape/rate hyperparameters for \eqn{\tau_k}
#'   (defaults 1 and 1/50).
#' @param a.sig,b.sig Inverse-gamma hyperparameters for the idiosyncratic
#'   variance \eqn{\sigma^2} (defaults 3, 3).
#' @param a.vs2,b.vs2 Inverse-gamma hyperparameters for the subject random-effect
#'   variance \eqn{u_s^2} (defaults 3, 3).
#' @param acc.target Target acceptance rate for the adaptive Metropolis update of
#'   \eqn{\phi_j} (default 0.44).
#' @param Lr,a.phi.r,a.w,b.w,ur2 Hyperparameters of the mean-constrained
#'   mixture-of-mixtures prior on the sample-size factors \eqn{r_i}: number of
#'   mixture components, Dirichlet concentration, two Beta parameters and the
#'   component variance.
#' @param L.theta,a.phi.theta,a.w.theta,b.w.theta,u2.theta Corresponding
#'   hyperparameters of the mixture-of-mixtures prior on the OTU-size factors
#'   \eqn{\alpha_j}.
#' @param m_kappa,sig_kappa Prior mean and variance for the probit
#'   zero-inflation coefficients \eqn{\kappa} (defaults 0, 3).
#' @param m_beta,sig_beta Prior mean and variance for the abundance regression
#'   coefficients \eqn{\beta} (defaults 0, 25). Ignored when `X = NULL`.
#'
#' @return A list of length `niter - burnin`. Each element is one posterior draw,
#'   itself a list with components:
#'   \describe{
#'     \item{`Lambda`}{`J` by `K` factor loading matrix.}
#'     \item{`sig2`}{idiosyncratic variance \eqn{\sigma^2}.}
#'     \item{`vs2`}{subject random-effect variance \eqn{u_s^2}.}
#'     \item{`phi`,`tau.k`}{Dirichlet-Laplace parameters.}
#'     \item{`eta`}{`n` by `K` latent factor scores.}
#'     \item{`ri`,`thetaj`}{sample-size factors \eqn{r_i} and OTU-size factors
#'       \eqn{\alpha_j} (only \eqn{r_i + \alpha_j} is identifiable).}
#'     \item{`sij`,`smj`}{sample- and subject-level random effects.}
#'     \item{`beta`}{`J` by `P` abundance regression coefficients (present only
#'       when `X` is supplied).}
#'     \item{`delta`}{`n` by `J` presence indicators; `1` = OTU present (drawn
#'       from the count distribution), `0` = structural/absent zero. (Note this
#'       is the complement of the \eqn{\delta_{ij}} in the paper.)}
#'     \item{`eps.ij`}{`n` by `J` matrix of estimated absence probabilities
#'       \eqn{\epsilon_{ij}} (the paper's zero-inflation probability).}
#'     \item{`kappa`}{`J` by `P + 1` probit zero-inflation coefficients
#'       (intercept in the first column).}
#'   }
#'
#' Common posterior summaries: the marginal OTU correlation matrix is
#' `cov2cor(tcrossprod(draw$Lambda) + diag(draw$vs2 + draw$sig2, J))`, averaged
#' over draws.
#'
#' @references
#' Zhang, S., Shuler, K., Lee, J. and Jeganathan, P. (2023). Bayesian Modeling of
#' Interaction between Features in Sparse Multivariate Count Data with Application
#' to Microbiome Study. *Annals of Applied Statistics*, 17(3).
#' \doi{10.1214/22-AOAS1690}
#'
#' @examples
#' \donttest{
#' sim <- simulate_zimln(n = 30, J = 30, K = 3, seed = 1)
#' fit <- ZI_MLN(sim$Y, m = sim$m, M = sim$M, niter = 400, burnin = 200)
#' # posterior mean marginal correlation
#' J <- ncol(sim$Y)
#' rho <- Reduce(`+`, lapply(fit, function(d)
#'   cov2cor(tcrossprod(d$Lambda) + diag(d$vs2 + d$sig2, J)))) / length(fit)
#' }
#'
#' @export
ZI_MLN <- function(Y, X = NULL, m = NULL, M = NULL, K = 10,
                   niter = 20000, burnin = floor(niter / 2), seed = 3,
                   a.phi = 1/2, a.tau = 1, b.tau = 1/50,
                   a.sig = 3, b.sig = 3, a.vs2 = 3, b.vs2 = 3,
                   acc.target = 0.44,
                   Lr = 8, a.phi.r = 3, a.w = 5, b.w = 5, ur2 = 1,
                   L.theta = 15, a.phi.theta = 1, a.w.theta = 5,
                   b.w.theta = 5, u2.theta = 1,
                   m_kappa = 0, sig_kappa = 3, m_beta = 0, sig_beta = 25) {

  Y <- as.matrix(Y)
  n <- nrow(Y)
  J <- ncol(Y)
  if (is.null(m)) m <- seq_len(n)
  if (length(m) != n) stop("`m` must have length nrow(Y).")
  if (is.null(M)) M <- length(unique(m))
  stopifnot(niter > burnin, burnin >= 0)

  has.X <- !is.null(X)
  if (has.X) {
    X <- as.matrix(X)
    if (nrow(X) != n) stop("`X` must have nrow(Y) rows.")
    p <- ncol(X)
  } else {
    p <- 0L
  }
  ## n x (p + 1) design for the probit; intercept-only when X is NULL
  ## (cbind(1, NULL) would collapse to 1x1, so build it explicitly)
  tilde.X <- if (has.X) cbind(1, X) else matrix(1, n, 1)

  set.seed(seed)

  ## --- initial values ---------------------------------------------------
  zeta.m <- matrix(rexp(J * K, 1/2), J, K)
  t.phi <- rgamma(J, 1, 1)
  t.t.phi <- log(t.phi)
  phi <- t.phi / sum(t.phi)
  tau.k <- rgamma(K, a.tau, b.tau)
  sig2 <- rgamma(1, a.sig, b.sig)
  vs2 <- 1 / rgamma(1, a.vs2, b.vs2)
  eta <- rmvnorm(n, rep(0, K), diag(1, K))
  Lambda <- matrix(0, J, K)
  sig.pro <- rep(sqrt(2.38^2 / J), J)
  nacc <- rep(0, J)
  acc <- rep(0, J)

  ## --- empirical mean-constraint targets v_r and v_alpha ----------------
  ql <- seq(0.5, 0.99, by = 0.01)
  ql.matrix <- sapply(1:n, function(x) quantile(Y[x, ], ql))
  dl.new <- vapply(seq_len(nrow(ql.matrix)), function(i) {
    ql.medj <- median(ql.matrix[i, ])
    median(abs(ql.matrix[i, ] - ql.medj))
  }, numeric(1))
  find.best.l <- which(dl.new[-1] > 1.1 * dl.new[-length(dl.new)])
  best.l <- ql[find.best.l[1]]
  hat.ri <- vapply(1:n, function(i)
    log(0.01 + sum(Y[i, Y[i, ] < quantile(Y[i, ], best.l)])), numeric(1))
  hat.thetaj <- (log(colSums(Y + 0.01)) - sum(hat.ri)) / n

  vr <- mean(hat.ri); a.xi <- vr; sig2.xi <- 1
  v.theta <- mean(hat.thetaj); a.xi.theta <- v.theta; sig2.xi.theta <- 1

  ## --- starting points for r_i / alpha_j mixtures -----------------------
  ri <- hat.ri
  Si1 <- sample(1:Lr, n, replace = TRUE)
  Si2 <- sample(0:1, n, replace = TRUE)
  xi <- rep(0, Lr)
  w.l <- rbeta(Lr, a.w, b.w)

  thetaj <- hat.thetaj
  Sj1 <- sample(1:L.theta, J, replace = TRUE)
  Sj2 <- sample(0:1, J, replace = TRUE)
  xi.theta <- rep(0, L.theta)
  w.l.theta <- rbeta(L.theta, a.w.theta, b.w.theta)

  y.star <- matrix(log(Y + 0.01), n, J)
  smj <- matrix(0, M, J)
  sij <- smj[m, ]

  ## --- zero-inflation (probit) initial values ---------------------------
  eps.ij <- matrix(runif(n * J), n, J)
  delta <- matrix(1, n, J)
  Z <- matrix(0, n, J)
  kappa <- matrix(0, J, p + 1)
  mu_kappa <- matrix(rep(m_kappa, p + 1), ncol = 1)
  Sigma_kappa_inv <- diag(1 / sig_kappa, p + 1)

  ## --- abundance regression initial values ------------------------------
  if (has.X) {
    beta <- matrix(0, J, p)
    mu_beta <- matrix(rep(m_beta, p), ncol = 1)
    Sigma_beta_inv <- diag(1 / sig_beta, p)
  }

  ## log target for the M-H update of phi_j (j passed explicitly)
  posterior.phi <- function(t.phi, j) {
    sum(dgamma(t.phi, a.phi, 1, log = TRUE)) +
      sum(dlaplace(c(Lambda), mu = 0,
                   sigma = rep(t.phi / sum(t.phi), K) * rep(tau.k, each = J),
                   log = TRUE)) +
      log(t.phi[j])
  }

  index.Y.0 <- which(Y == 0, arr.ind = TRUE)
  index.Y.non0 <- which(Y > 0, arr.ind = TRUE)
  delta[index.Y.non0] <- 1
  delta[index.Y.0] <- 0

  RSS <- y.star - matrix(ri, n, J) - matrix(thetaj, n, J, byrow = TRUE) -
    sij - tcrossprod(eta, Lambda)
  if (has.X) RSS <- RSS - tcrossprod(X, beta)
  logY <- log(Y)
  logY.1 <- log(Y + 1)

  ls <- vector("list", niter - burnin)

  for (i in 1:niter) {

    ## ----- update delta (presence indicator) for observed zeros ---------
    RSS <- y.star - RSS                      # RSS now holds the fitted mean
    ## vapply (not sapply): with zero matching entries sapply returns list(),
    ## which would coerce the matrix to a list and destroy its dimensions.
    delta[index.Y.0] <- vapply(seq_len(nrow(index.Y.0)), function(x) {
      prob.0 <- eps.ij[index.Y.0[x, 1], index.Y.0[x, 2]]
      prob.1 <- (1 - prob.0) *
        pnorm(0, mean = RSS[index.Y.0[x, 1], index.Y.0[x, 2]], sd = sqrt(sig2))
      as.numeric(rbinom(1, 1, prob = prob.1 / (prob.0 + prob.1)))
    }, numeric(1))
    index.delta.1 <- which(delta == 1, arr.ind = TRUE)   # present
    index.delta.0 <- which(delta == 0, arr.ind = TRUE)   # absent
    RSS <- y.star - RSS                      # back to residual

    ## ----- probit latent Z ----------------------------------------------
    Z[index.delta.1] <- vapply(seq_len(nrow(index.delta.1)), function(x)
      rtruncnorm(1, a = -Inf, b = 0,
                 mean = tilde.X[index.delta.1[x, 1], ] %*% kappa[index.delta.1[x, 2], ],
                 sd = 1), numeric(1))
    Z[index.delta.0] <- vapply(seq_len(nrow(index.delta.0)), function(x)
      rtruncnorm(1, a = 0, b = Inf,
                 mean = tilde.X[index.delta.0[x, 1], ] %*% kappa[index.delta.0[x, 2], ],
                 sd = 1), numeric(1))

    ## ----- update kappa and epsilon -------------------------------------
    pos.var.kappa <- solve(crossprod(tilde.X) + Sigma_kappa_inv)
    kappa <- t(sapply(1:J, function(x)
      rmvn(1, pos.var.kappa %*% (Sigma_kappa_inv %*% mu_kappa +
                                   crossprod(tilde.X, Z[, x, drop = FALSE])),
           pos.var.kappa)))
    if (p == 0) kappa <- matrix(kappa, J, 1)
    eps.ij <- pnorm(tcrossprod(tilde.X, kappa))

    ## ----- subject random effects s_mj ----------------------------------
    RSS <- RSS + sij
    for (msample in 1:M) {
      for (j in 1:J) {
        isample <- which(m == msample)
        if (length(which(delta[isample, j] == 1)) == 0) {
          smj[msample, j] <- rnorm(1) * sqrt(vs2)
        } else {
          isample <- isample[which(delta[isample, j] == 1)]
          smj.var <- 1 / (1 / vs2 + length(isample) / sig2)
          smj.mean <- smj.var * sum(RSS[isample, j]) / sig2
          smj[msample, j] <- rnorm(1) * sqrt(smj.var) + smj.mean
        }
      }
    }
    sij <- smj[m, ]
    RSS <- RSS - sij

    ## ----- update u_s^2 (vs2) -------------------------------------------
    vs2 <- 1 / rgamma(1, a.vs2 + M * J / 2, b.vs2 + sum(smj^2) / 2)

    ## ----- impute continuous latent y* on present entries ---------------
    RSS <- y.star - RSS
    mean_uij <- RSS[index.delta.1]
    u_min <- pnorm(logY[index.delta.1], mean = mean_uij, sd = sqrt(sig2))
    u_max <- pnorm(logY.1[index.delta.1], mean = mean_uij, sd = sqrt(sig2))
    sample_u <- runif(nrow(index.delta.1), u_min, u_max)
    sample_u[sample_u == 1] <- 1 - 9.9e-16
    sample_u[sample_u == 0] <- 9.9e-16
    y.star[index.delta.1] <- qnorm(sample_u, mean = mean_uij, sd = sqrt(sig2))
    RSS <- y.star - RSS

    ## ----- Lambda -------------------------------------------------------
    RSS <- RSS + tcrossprod(eta, Lambda)
    for (j in 1:J) {
      index.lambda <- which(delta[, j] != 0)
      Q <- diag(1 / zeta.m[j, ] * 1 / phi[j]^2 * 1 / tau.k^2) +
        crossprod(eta[index.lambda, , drop = FALSE]) / sig2
      U <- chol(Q)
      a <- crossprod(eta[index.lambda, , drop = FALSE], RSS[index.lambda, j]) / sig2
      Lambda[j, ] <- backsolve(U, backsolve(U, a, transpose = TRUE) + rnorm(K))
    }

    ## ----- eta ----------------------------------------------------------
    for (isample in 1:n) {
      index.etai <- which(delta[isample, ] != 0)
      ## drop = FALSE: a sample with a single present OTU would otherwise
      ## collapse Lambda[index.etai, ] to a vector and break crossprod()
      Q <- diag(1, nrow = K) + crossprod(Lambda[index.etai, , drop = FALSE]) / sig2
      a <- crossprod(Lambda[index.etai, , drop = FALSE], RSS[isample, index.etai]) / sig2
      U <- chol(Q)
      eta[isample, ] <- backsolve(U, backsolve(U, a, transpose = TRUE) + rnorm(K))
    }
    RSS <- RSS - tcrossprod(eta, Lambda)

    ## ----- sigma^2 ------------------------------------------------------
    inter.RSS <- RSS
    inter.RSS[which(delta == 0, arr.ind = TRUE)] <- 0
    sig2 <- 1 / rgamma(1, a.sig + sum(delta) / 2,
                       rate = b.sig + 1/2 * sum(inter.RSS^2))

    ## ----- Dirichlet-Laplace augmentation: zeta, tau, phi ---------------
    mean.psi <- matrix(phi, J, K) * matrix(tau.k, J, K, byrow = TRUE) / abs(Lambda)
    zeta.m <- matrix(1 / rinvgauss(J * K, mean = c(mean.psi), dispersion = 1), J, K)

    chi.vec <- colSums(2 * abs(Lambda) / matrix(phi, J, K))
    tau.k <- sapply(1:K, function(x)
      rgig(1, lambda = a.tau - J, psi = 2 * b.tau, chi = chi.vec[x]))

    d <- min(0.01, 1 / sqrt(i))
    for (j in 1:J) {
      new.pro <- rnorm(1) * sig.pro[j] + t.t.phi[j]
      t.phi.pro <- t.phi
      t.phi.pro[j] <- exp(new.pro)
      d.log <- posterior.phi(t.phi.pro, j) - posterior.phi(t.phi, j)
      acc[j] <- min(1, exp(d.log))
      if (runif(1) < acc[j]) {
        t.t.phi[j] <- new.pro
        t.phi[j] <- exp(new.pro)
        phi <- t.phi / sum(t.phi)
        nacc[j] <- nacc[j] + 1
      }
      if (i %% 50 == 0) {
        ind <- ifelse(nacc[j] / 50 > acc.target, 1, -1)
        sig.pro[j] <- sig.pro[j] * exp(d * ind)
        nacc[j] <- 0
      }
    }

    ## ----- sample-size factors r_i --------------------------------------
    RSS <- RSS + matrix(ri, n, J)
    pos.v.ri <- 1 / (1 / ur2 + rowSums(delta) / sig2)
    prior.m.ri <- xi[Si1] * (Si2 == 1) +
      ((vr - w.l[Si1] * xi[Si1]) / (1 - w.l[Si1])) * (Si2 == 0)
    sample.mean <- sapply(1:n, function(x) sum(RSS[x, which(delta[x, ] != 0)]))
    ri <- rnorm(n) * sqrt(pos.v.ri) + (prior.m.ri / ur2 + sample.mean / sig2) * pos.v.ri
    RSS <- RSS - matrix(ri, n, J)

    phi.r <- rdirichlet(1, sapply(1:Lr, function(x) sum(Si1 == x)) + a.phi.r)
    w.l <- sapply(1:Lr, function(l)
      rbeta(1, sum(Si1 == l & Si2 == 1) + a.w, sum(Si1 == l & Si2 == 0) + b.w))
    xi.pos.var <- sapply(1:Lr, function(l)
      1 / (1 / sig2.xi + sum(Si1 == l & Si2 == 1) / ur2 +
             w.l[l]^2 * sum(Si1 == l & Si2 == 0) / (ur2 * (1 - w.l[l])^2)))
    xi.pos.mean <- sapply(1:Lr, function(l)
      xi.pos.var[l] * (a.xi / sig2.xi + sum(ri[Si1 == l & Si2 == 1]) / ur2 +
        w.l[l] * (sum(Si1 == l & Si2 == 0) * vr -
                    (1 - w.l[l]) * sum(ri[Si1 == l & Si2 == 0])) /
          (ur2 * (1 - w.l[l])^2)))
    xi <- rnorm(Lr) * sqrt(xi.pos.var) + xi.pos.mean

    for (isample in 1:n) {
      pil1 <- dnorm(ri[isample], xi, sqrt(ur2)) * w.l * phi.r
      pil0 <- dnorm(ri[isample], (vr - w.l * xi) / (1 - w.l), sqrt(ur2)) *
        (1 - w.l) * phi.r
      pi.normalized <- c(pil1, pil0) / sum(pil1 + pil0)
      ind <- sample(1:(2 * Lr), 1, prob = pi.normalized)
      Si2[isample] <- ifelse(ind > Lr, 0, 1)
      Si1[isample] <- ifelse(ind > Lr, ind - Lr, ind)
    }

    ## ----- OTU-size factors alpha_j (thetaj) ----------------------------
    RSS <- RSS + matrix(thetaj, n, J, byrow = TRUE)
    pos.v.theta <- 1 / (1 / u2.theta + colSums(delta) / sig2)
    thetaj.prior.mean <- xi.theta[Sj1] * (Sj2 == 1) +
      ((v.theta - w.l.theta[Sj1] * xi.theta[Sj1]) / (1 - w.l.theta[Sj1])) * (Sj2 == 0)
    theta.sample.mean <- sapply(1:J, function(x) sum(RSS[which(delta[, x] != 0), x]))
    thetaj <- rnorm(J) * sqrt(pos.v.theta) +
      (thetaj.prior.mean / u2.theta + theta.sample.mean / sig2) * pos.v.theta
    RSS <- RSS - matrix(thetaj, n, J, byrow = TRUE)

    phi.r.theta <- rdirichlet(1, sapply(1:L.theta, function(x) sum(Sj1 == x)) + a.phi.theta)
    w.l.theta <- sapply(1:L.theta, function(l)
      rbeta(1, sum(Sj1 == l & Sj2 == 1) + a.w.theta,
            sum(Sj1 == l & Sj2 == 0) + b.w.theta))
    xi.pos.var <- sapply(1:L.theta, function(l)
      1 / (1 / sig2.xi.theta + sum(Sj1 == l & Sj2 == 1) / u2.theta +
             w.l.theta[l]^2 * sum(Sj1 == l & Sj2 == 0) / (u2.theta * (1 - w.l.theta[l])^2)))
    xi.pos.mean <- sapply(1:L.theta, function(l)
      xi.pos.var[l] * (a.xi.theta / sig2.xi.theta +
        sum(thetaj[Sj1 == l & Sj2 == 1]) / u2.theta +
        w.l.theta[l] * (sum(Sj1 == l & Sj2 == 0) * v.theta -
                          (1 - w.l.theta[l]) * sum(thetaj[Sj1 == l & Sj2 == 0])) /
          (u2.theta * (1 - w.l.theta[l])^2)))
    xi.theta <- rnorm(L.theta, xi.pos.mean, sqrt(xi.pos.var))

    for (jsample in 1:J) {
      pil1 <- dnorm(thetaj[jsample], xi.theta, sqrt(u2.theta)) * w.l.theta * phi.r.theta
      pil0 <- dnorm(thetaj[jsample], (v.theta - w.l.theta * xi.theta) / (1 - w.l.theta),
                    sqrt(u2.theta)) * (1 - w.l.theta) * phi.r.theta
      pi.normalized <- c(pil1, pil0) / sum(pil1 + pil0)
      ind <- sample(1:(2 * L.theta), 1, prob = pi.normalized)
      Sj1[jsample] <- ifelse(ind > L.theta, ind - L.theta, ind)
      Sj2[jsample] <- ifelse(ind > L.theta, 0, 1)
    }

    ## ----- abundance regression beta ------------------------------------
    if (has.X) {
      RSS <- RSS + tcrossprod(X, beta)
      for (j in 1:J) {
        index.betaj <- which(delta[, j] != 0)
        Q <- Sigma_beta_inv + crossprod(X[index.betaj, , drop = FALSE]) / sig2
        U <- chol(Q)
        a <- Sigma_beta_inv %*% mu_beta +
          crossprod(X[index.betaj, , drop = FALSE], RSS[index.betaj, j]) / sig2
        beta[j, ] <- backsolve(U, backsolve(U, a, transpose = TRUE) + rnorm(p))
      }
      RSS <- RSS - tcrossprod(X, beta)
    }

    ## ----- store draw ---------------------------------------------------
    if (i > burnin) {
      draw <- list(Lambda = Lambda, sig2 = sig2, vs2 = vs2, phi = phi,
                   tau.k = tau.k, eta = eta, ri = ri, thetaj = thetaj,
                   sij = sij, smj = smj, delta = delta, eps.ij = eps.ij,
                   kappa = kappa)
      if (has.X) draw$beta <- beta
      ls[[i - burnin]] <- draw
    }
  }
  ls
}
