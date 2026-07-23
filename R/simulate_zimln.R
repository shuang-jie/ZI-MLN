#' Simulate a count table from the ZI-MLN generative model
#'
#' Draws a synthetic OTU count matrix (and, optionally, covariates) directly from
#' the zero-inflated multivariate rounded log-normal model, together with the
#' true parameter values used to generate it. Useful for examples, unit tests and
#' checking parameter recovery with [ZI_MLN()].
#'
#' @param n Number of samples.
#' @param J Number of OTUs (features).
#' @param K Number of latent factors used to build the true covariance.
#' @param zero.rate Proportion of entries in each column of the true loading
#'   matrix `Lambda` set to zero (controls covariance sparsity). Default 0.8.
#' @param sig2 True idiosyncratic variance \eqn{\sigma^2} (default 1).
#' @param vs2 True subject random-effect variance \eqn{u_s^2} (default 1).
#' @param p Number of continuous covariates to generate. `0` (default) produces a
#'   no-covariate data set (`X` is `NULL`).
#' @param M Number of subjects/groups (default `n`, i.e. one sample per subject).
#' @param m Optional length-`n` subject index. Defaults to a balanced assignment
#'   of the `n` samples to the `M` subjects.
#' @param seed Random seed (default 1).
#'
#' @return A list with the simulated data and the ground truth:
#'   \describe{
#'     \item{`Y`}{`n` by `J` integer count matrix (the model input).}
#'     \item{`X`}{`n` by `p` covariate matrix, or `NULL` when `p = 0`.}
#'     \item{`m`,`M`}{subject index and subject count.}
#'     \item{`Lambda`,`sig2`,`vs2`}{true covariance components.}
#'     \item{`Omega`,`true.cor`}{true covariance `Lambda Lambda' + sig2 I` and the
#'       true marginal correlation matrix `cov2cor(Omega + vs2 I)`.}
#'     \item{`ri`,`thetaj`,`smj`}{true sample-size / OTU-size factors and subject
#'       random effects.}
#'     \item{`beta`}{true `J` by `p` regression coefficients (`NULL` when `p = 0`).}
#'     \item{`kappa`}{true `J` by `p + 1` probit zero-inflation coefficients.}
#'     \item{`eps`}{`n` by `J` true absence probabilities \eqn{\epsilon_{ij}}.}
#'     \item{`delta`}{`n` by `J` true presence indicators (`1` = present).}
#'   }
#'
#' @seealso [ZI_MLN()]
#'
#' @examples
#' sim <- simulate_zimln(n = 20, J = 30, K = 3, seed = 1)
#' dim(sim$Y)
#' mean(sim$Y == 0)   # overall zero rate
#'
#' @export
simulate_zimln <- function(n = 20, J = 150, K = 5, zero.rate = 0.8,
                           sig2 = 1, vs2 = 1, p = 0, M = n, m = NULL,
                           seed = 1) {
  set.seed(seed)
  if (is.null(m)) m <- rep(seq_len(M), length.out = n)
  if (length(m) != n) stop("`m` must have length n.")

  ## true loading matrix with column-wise sparsity
  Lambda <- matrix(runif(J * K, -3, 3), J, K)
  for (k in 1:K) {
    ind <- sample(1:J, size = round(J * zero.rate))
    Lambda[ind, k] <- 0
  }
  Omega <- tcrossprod(Lambda) + sig2 * diag(J)

  ri <- runif(n, 3, 7)
  thetaj <- runif(J, 0, 2)
  smj <- matrix(rnorm(M * J, 0, sqrt(vs2)), M, J)

  ## covariates and abundance regression
  if (p > 0) {
    X <- matrix(rnorm(n * p), n, p)
    beta <- matrix(rnorm(J * p, 0, 1), J, p)
    Xbeta <- tcrossprod(X, beta)
  } else {
    X <- NULL
    beta <- NULL
    Xbeta <- matrix(0, n, J)
  }
  tilde.X <- if (p > 0) cbind(1, X) else matrix(1, n, 1)

  mu <- matrix(ri, n, J) + matrix(thetaj, n, J, byrow = TRUE) +
    smj[m, , drop = FALSE] + Xbeta

  y.star <- t(sapply(1:n, function(i) rmvnorm(1, mean = mu[i, ], sigma = Omega)))
  Y <- floor(exp(y.star))

  ## zero inflation (probit): eps = P(absent)
  kappa <- matrix(runif(J * (p + 1), -1, 0), J, p + 1)
  eps <- pnorm(tcrossprod(tilde.X, kappa))
  delta <- matrix(rbinom(n * J, 1, 1 - eps), n, J)   # 1 = present
  Y <- Y * delta

  list(Y = Y, X = X, m = m, M = M,
       Lambda = Lambda, sig2 = sig2, vs2 = vs2,
       Omega = Omega, true.cor = cov2cor(Omega + diag(vs2, J)),
       ri = ri, thetaj = thetaj, smj = smj, beta = beta,
       kappa = kappa, eps = eps, delta = delta)
}
