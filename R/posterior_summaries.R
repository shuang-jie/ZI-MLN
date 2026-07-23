## Posterior summaries of a ZI_MLN() fit (a list of MCMC draws).

#' Posterior marginal correlation between OTUs
#'
#' Summarises the posterior of the marginal correlation matrix
#' \eqn{\rho_{jj'}} implied by a [ZI_MLN()] fit. For each retained draw the
#' marginal correlation is computed from the covariance
#' \eqn{\Omega = \Lambda\Lambda' + (\sigma^2 + u_s^2) I} as
#' `cov2cor(tcrossprod(Lambda) + diag(vs2 + sig2, J))`, and the draws are then
#' summarised element-wise.
#'
#' @param fit A list returned by [ZI_MLN()] (one element per posterior draw).
#' @param ci Logical; if `TRUE` also return element-wise credible bounds. This
#'   builds a `J x J x ndraw` array and can use substantial memory for large `J`.
#'   Default `FALSE` (posterior mean only).
#' @param prob Length-2 vector of lower/upper probabilities for the credible
#'   interval when `ci = TRUE`. Default `c(0.025, 0.975)`.
#'
#' @return If `ci = FALSE`, the `J x J` posterior-mean correlation matrix. If
#'   `ci = TRUE`, a list with `mean`, `lower` and `upper` `J x J` matrices.
#'
#' @seealso [posterior_covariance()], [posterior_beta()]
#' @examples
#' \donttest{
#' sim <- simulate_zimln(n = 20, J = 20, K = 2, seed = 1)
#' fit <- ZI_MLN(sim$Y, m = sim$m, M = sim$M, niter = 400, burnin = 200)
#' rho <- posterior_correlation(fit)
#' }
#' @export
posterior_correlation <- function(fit, ci = FALSE, prob = c(0.025, 0.975)) {
  .check_fit(fit)
  J <- nrow(fit[[1]]$Lambda)
  draw_mat <- function(d) cov2cor(tcrossprod(d$Lambda) + diag(d$vs2 + d$sig2, J))
  .summarise_matrix(fit, draw_mat, ci, prob)
}

#' Posterior covariance between OTUs
#'
#' Summarises the posterior of the OTU covariance matrix from a [ZI_MLN()] fit.
#' By default the *marginal* covariance
#' \eqn{\Omega = \Lambda\Lambda' + (\sigma^2 + u_s^2) I} is returned (the one
#' whose correlation is [posterior_correlation()]). Set `marginal = FALSE` for the
#' interaction covariance \eqn{\Sigma = \Lambda\Lambda' + \sigma^2 I}, which
#' excludes the subject random-effect variance.
#'
#' @param fit A list returned by [ZI_MLN()].
#' @param marginal Logical; add the subject random-effect variance \eqn{u_s^2} to
#'   the diagonal (`TRUE`, default, giving \eqn{\Omega}) or not (`FALSE`, giving
#'   \eqn{\Sigma}).
#' @param ci,prob As in [posterior_correlation()].
#'
#' @return A `J x J` posterior-mean covariance matrix, or a list of
#'   `mean`/`lower`/`upper` matrices when `ci = TRUE`.
#'
#' @seealso [posterior_correlation()], [posterior_beta()]
#' @examples
#' \donttest{
#' sim <- simulate_zimln(n = 20, J = 20, K = 2, seed = 1)
#' fit <- ZI_MLN(sim$Y, m = sim$m, M = sim$M, niter = 400, burnin = 200)
#' Sigma <- posterior_covariance(fit, marginal = FALSE)
#' }
#' @export
posterior_covariance <- function(fit, marginal = TRUE, ci = FALSE,
                                 prob = c(0.025, 0.975)) {
  .check_fit(fit)
  J <- nrow(fit[[1]]$Lambda)
  draw_mat <- function(d) {
    add <- if (marginal) d$vs2 + d$sig2 else d$sig2
    tcrossprod(d$Lambda) + diag(add, J)
  }
  .summarise_matrix(fit, draw_mat, ci, prob)
}

#' Posterior summary of the abundance regression coefficients beta
#'
#' Summarises the posterior of the covariate effects \eqn{\beta_{jp}} from a
#' [ZI_MLN()] fit that was run *with* covariates. Errors if the model was fit
#' without covariates (`X = NULL`).
#'
#' @param fit A list returned by [ZI_MLN()] with a `beta` component in each draw.
#' @param prob Length-2 vector of lower/upper probabilities for the credible
#'   interval. Default `c(0.025, 0.975)`.
#'
#' @return A list with three `J x P` matrices: `mean`, `lower` and `upper`
#'   (posterior mean and credible bounds of each \eqn{\beta_{jp}}), plus a tidy
#'   long `data.frame` `table` with columns `otu`, `covariate`, `mean`, `lower`,
#'   `upper`.
#'
#' @seealso [posterior_correlation()], [posterior_covariance()]
#' @examples
#' \donttest{
#' sim <- simulate_zimln(n = 30, J = 20, K = 2, p = 2, seed = 1)
#' fit <- ZI_MLN(sim$Y, X = sim$X, m = sim$m, M = sim$M, niter = 400, burnin = 200)
#' b <- posterior_beta(fit)
#' head(b$table)
#' }
#' @export
posterior_beta <- function(fit, prob = c(0.025, 0.975)) {
  .check_fit(fit)
  if (is.null(fit[[1]]$beta))
    stop("This fit has no `beta`; the model was run without covariates (X = NULL).")
  J <- nrow(fit[[1]]$beta)
  P <- ncol(fit[[1]]$beta)
  arr <- vapply(fit, function(d) d$beta, matrix(0, J, P))  # J x P x ndraw
  est   <- apply(arr, 1:2, mean)
  lower <- apply(arr, 1:2, quantile, probs = prob[1])
  upper <- apply(arr, 1:2, quantile, probs = prob[2])
  dimnames(est) <- dimnames(lower) <- dimnames(upper) <-
    list(otu = seq_len(J), covariate = seq_len(P))
  table <- data.frame(
    otu = rep(seq_len(J), times = P),
    covariate = rep(seq_len(P), each = J),
    mean = c(est), lower = c(lower), upper = c(upper)
  )
  list(mean = est, lower = lower, upper = upper, table = table)
}

## ---- internal helpers -----------------------------------------------------

.check_fit <- function(fit) {
  if (!is.list(fit) || length(fit) == 0 || is.null(fit[[1]]$Lambda))
    stop("`fit` must be a non-empty list of draws returned by ZI_MLN().")
  invisible(TRUE)
}

## Posterior mean (and optional element-wise CI) of a per-draw matrix.
.summarise_matrix <- function(fit, draw_mat, ci, prob) {
  B <- length(fit)
  if (!ci) {
    m <- draw_mat(fit[[1]])
    for (b in seq_len(B)[-1]) m <- m + draw_mat(fit[[b]])  # incremental: no big array
    return(m / B)
  }
  arr <- vapply(fit, draw_mat, draw_mat(fit[[1]]))  # J x J x B
  list(mean  = apply(arr, 1:2, mean),
       lower = apply(arr, 1:2, quantile, probs = prob[1]),
       upper = apply(arr, 1:2, quantile, probs = prob[2]))
}
