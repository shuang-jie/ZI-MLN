## Fast structural tests (run on every R CMD check).

test_that("simulate_zimln returns coherent shapes", {
  s0 <- simulate_zimln(n = 12, J = 10, K = 2, p = 0, seed = 1)
  expect_equal(dim(s0$Y), c(12, 10))
  expect_null(s0$X)
  expect_null(s0$beta)
  expect_true(all(s0$delta %in% c(0, 1)))
  expect_true(all(s0$Y[s0$delta == 0] == 0))          # absent => zero count

  s1 <- simulate_zimln(n = 12, J = 10, K = 2, p = 2, seed = 1)
  expect_equal(dim(s1$X), c(12, 2))
  expect_equal(dim(s1$beta), c(10, 2))
  expect_equal(dim(s1$kappa), c(10, 3))               # p + 1 columns
})

test_that("ZI_MLN (no covariate) runs and returns well-formed draws", {
  s <- simulate_zimln(n = 15, J = 12, K = 2, p = 0, seed = 2)
  fit <- ZI_MLN(s$Y, m = s$m, M = s$M, K = 4, niter = 40, burnin = 20, seed = 2)

  expect_length(fit, 20)
  d <- fit[[1]]
  expect_equal(dim(d$Lambda), c(12, 4))
  expect_equal(dim(d$eta), c(15, 4))
  expect_equal(dim(d$eps.ij), c(15, 12))
  expect_equal(dim(d$delta), c(15, 12))
  expect_equal(dim(d$kappa), c(12, 1))                # intercept only
  expect_null(d$beta)

  expect_true(all(d$eps.ij > 0 & d$eps.ij < 1))
  expect_true(all(d$delta %in% c(0, 1)))
  expect_true(all(d$delta[s$Y > 0] == 1))             # positive counts present
  expect_false(anyNA(d$Lambda))
})

test_that("ZI_MLN (with covariate) runs and returns beta and kappa columns", {
  s <- simulate_zimln(n = 16, J = 12, K = 2, p = 2, M = 8,
                      m = rep(1:8, 2), seed = 3)
  fit <- ZI_MLN(s$Y, X = s$X, m = s$m, M = s$M, K = 4,
                niter = 40, burnin = 20, seed = 3)

  expect_length(fit, 20)
  d <- fit[[1]]
  expect_equal(dim(d$beta), c(12, 2))
  expect_equal(dim(d$kappa), c(12, 3))                # intercept + 2 covariates
  expect_false(anyNA(d$beta))
})

test_that("posterior summary functions return correct shapes", {
  s <- simulate_zimln(n = 16, J = 10, K = 2, p = 2, M = 8,
                      m = rep(1:8, 2), seed = 5)
  fit <- ZI_MLN(s$Y, X = s$X, m = s$m, M = s$M, K = 3,
                niter = 40, burnin = 20, seed = 5)

  rho <- posterior_correlation(fit)
  expect_equal(dim(rho), c(10, 10))
  expect_true(all(abs(rho) <= 1 + 1e-8))
  expect_equal(diag(rho), rep(1, 10), tolerance = 1e-6)

  rho.ci <- posterior_correlation(fit, ci = TRUE)
  expect_named(rho.ci, c("mean", "lower", "upper"))
  expect_true(all(rho.ci$lower <= rho.ci$upper + 1e-8))

  Sig <- posterior_covariance(fit, marginal = FALSE)
  Om  <- posterior_covariance(fit, marginal = TRUE)
  expect_equal(dim(Sig), c(10, 10))
  expect_true(all(diag(Om) >= diag(Sig) - 1e-8))   # marginal adds vs2 to diagonal

  b <- posterior_beta(fit)
  expect_equal(dim(b$mean), c(10, 2))
  expect_equal(nrow(b$table), 20)
  expect_true(all(b$table$lower <= b$table$upper + 1e-8))
})

test_that("posterior_beta errors without covariates", {
  s <- simulate_zimln(n = 12, J = 8, K = 2, p = 0, seed = 6)
  fit <- ZI_MLN(s$Y, m = s$m, M = s$M, K = 3, niter = 30, burnin = 15, seed = 6)
  expect_error(posterior_beta(fit), "without covariates")
})

test_that("a near-empty sample (single present OTU) does not error", {
  # Regression test: Lambda[index.etai, ] used to drop to a vector when a sample
  # had very few present OTUs, breaking crossprod() in the eta update. The rest
  # of the table stays informative so the sampler itself remains well behaved.
  s <- simulate_zimln(n = 14, J = 12, K = 2, seed = 4)
  Y <- s$Y
  Y[1, ] <- 0L; Y[1, 3] <- 7L        # sample 1 has a single non-zero count
  Y[2, ] <- 0L                       # sample 2 is entirely zero
  for (sd in c(1, 2, 3)) {
    expect_no_error(
      ZI_MLN(Y, m = s$m, M = s$M, K = 3, niter = 60, burnin = 30, seed = sd)
    )
  }
})

test_that("input validation errors are raised", {
  s <- simulate_zimln(n = 10, J = 8, K = 2, seed = 4)
  expect_error(ZI_MLN(s$Y, m = 1:3), "length")
  expect_error(ZI_MLN(s$Y, X = matrix(0, 3, 2)), "rows")
})
