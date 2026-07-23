## Full parameter-recovery test. Slow (~4 min), so skipped by default.
## Run it explicitly with:  ZIMLN_RUN_RECOVERY=1 Rscript -e 'testthat::test_local()'
## or via devtools::test(filter = "recovery") after setting the env var.

auc_fn <- function(lab, pred) {
  r <- rank(pred); n1 <- sum(lab == 1); n0 <- sum(lab == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[lab == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

post_cor <- function(fit, J) {
  Reduce(`+`, lapply(fit, function(d)
    cov2cor(tcrossprod(d$Lambda) + diag(d$vs2 + d$sig2, J)))) / length(fit)
}

test_that("no-covariate model recovers key parameters", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("ZIMLN_RUN_RECOVERY")),
              "set ZIMLN_RUN_RECOVERY=1 to run the slow recovery test")

  s <- simulate_zimln(n = 40, J = 40, K = 3, zero.rate = 0.5, seed = 101)
  fit <- ZI_MLN(s$Y, m = s$m, M = s$M, niter = 4000, seed = 7)
  J <- ncol(s$Y)

  rt.est <- Reduce(`+`, lapply(fit, function(d) outer(d$ri, d$thetaj, `+`))) / length(fit)
  expect_gt(cor(c(rt.est), c(outer(s$ri, s$thetaj, `+`))), 0.88)

  mc <- post_cor(fit, J); ut <- upper.tri(mc)
  expect_gt(cor(mc[ut], s$true.cor[ut]), 0.75)
  expect_lt(mean(abs(mc[ut] - s$true.cor[ut])), 0.15)

  tot <- mean(sapply(fit, function(d) d$sig2 + d$vs2))
  expect_lt(abs(tot - (s$sig2 + s$vs2)), 0.6)

  delta.pm <- Reduce(`+`, lapply(fit, function(d) d$delta)) / length(fit)
  zi <- which(s$Y == 0)
  expect_gt(auc_fn(s$delta[zi], delta.pm[zi]), 0.80)
})

test_that("covariate model recovers beta and structure", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("ZIMLN_RUN_RECOVERY")),
              "set ZIMLN_RUN_RECOVERY=1 to run the slow recovery test")

  s <- simulate_zimln(n = 60, J = 40, K = 3, zero.rate = 0.5, p = 2,
                      M = 30, m = rep(1:30, 2), seed = 202)
  fit <- ZI_MLN(s$Y, X = s$X, m = s$m, M = s$M, niter = 4000, seed = 7)
  J <- ncol(s$Y)

  mc <- post_cor(fit, J); ut <- upper.tri(mc)
  expect_gt(cor(mc[ut], s$true.cor[ut]), 0.80)

  beta.est <- Reduce(`+`, lapply(fit, function(d) d$beta)) / length(fit)
  expect_gt(cor(c(beta.est), c(s$beta)), 0.85)

  delta.pm <- Reduce(`+`, lapply(fit, function(d) d$delta)) / length(fit)
  zi <- which(s$Y == 0)
  expect_gt(auc_fn(s$delta[zi], delta.pm[zi]), 0.85)

  ## eps.ij must estimate P(ABSENT), not P(present)
  eps.est <- Reduce(`+`, lapply(fit, function(d) d$eps.ij)) / length(fit)
  expect_gt(cor(c(eps.est), c(s$eps)), 0.80)
})
