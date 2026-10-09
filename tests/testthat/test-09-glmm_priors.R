context("bglmer numerical results with fixef and cov priors")

source(system.file("common", "glmmData.R", package = "blme"))
control <- glmerControl(optimizer = "Nelder_Mead")
lme4Version <- packageVersion("lme4")

test_that("bglmer fits test data with fixef prior, matching previous version", {
  fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control,
                cov.prior = NULL, fixef.prior = normal)
  if (lme4Version < "1.1-4") {
    expect_equal(fit@theta, 1.26501253837861)
    expect_equal(fit@beta, c(0.873121247636467, 2.46647249930796, 1.32070156863358))
  } else if (lme4Version < "2.1-0") {
    expect_equal(fit@theta, 1.26501385482573)
    expect_equal(fit@beta, c(0.873131216784199, 2.46647899095567, 1.32070496549675))
  } else {
    expect_equal(fit@theta, 1.26501476305925, tolerance = 1.0e-6)
    expect_equal(fit@beta, c(0.873129799640665, 2.4664666989276, 1.32070290336581), tolerance = 1.0e-6)
  }
})

test_that("bglmer fits test data with cov prior, matching previous version", {
  fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control,
                cov.prior = wishart)
  if (lme4Version < "1.1-4") {
    expect_equal(fit@theta, 2.96766525351892)
    expect_equal(fit@beta, c(1.0963789854971, 3.67790570859986, 1.75655010020603))
  } else if (lme4Version < "2.1-0") {
    expect_equal(fit@theta, 2.96767284827046)
    expect_equal(fit@beta, c(1.0963789854971, 3.67790570859986, 1.75655010020603))
  } else {
    expect_equal(fit@theta, 2.96764095859915, tolerance = 1.0e-6)
    expect_equal(fit@beta, c(1.096412380265, 3.67792366613539, 1.75655008627344), tolerance = 1.0e-6)
  }
})

test_that("bglmer runs with fixef horseshoe prior", {
  skip_if_not_installed("expint")
  suppressWarnings(fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control,
                                 cov.prior = NULL, fixef.prior = horseshoe))
  
  ## the horseshoe density has a pole at zero, so the mode is sensitive to where
  ## the optimizer stops; check that the fit is sane rather than exact values
  expect_is(fit, "bglmerMod")
  expect_true(all(is.finite(c(fit@theta, fit@beta))))
  expect_lt(abs(fit@beta[1L]), 1e-4)
})

test_that("bglmer runs with fixef independent_ts prior", {
  fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control,
                cov.prior = NULL, fixef.prior = independent_ts(df = 1))
  
  expect_is(fit, "bglmerMod")
  expect_is(fit@priors$fixefPrior, "bmerIndependentTsDist")
  expect_true(all(is.finite(c(fit@theta, fit@beta))))
})

test_that("bglmer runs cov prior specified using level.dim", {
  fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control,
                cov.prior = wishart(df = level.dim + 1.25))
  expect_is(fit, "bglmerMod")
  expect_equal(fit@priors$cov[[1]]@df, 2.25)
})

test_that("blmer fits test data with custom prior, matching builtin gamma", {
  dgamma_cust <- function(x) dgamma(x, shape = 2.5, rate = 0.01, log = TRUE)
  
  fit.prof <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, control = control,
                     family = binomial(),
                     cov.prior = gamma(rate = 0.01))
  fit.cust <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, control = control,
                     family = binomial(),
                     cov.prior = custom(dgamma_cust, chol = TRUE, scale = "log"))
  expect_equal(fit.prof@theta, fit.cust@theta, tolerance = 1e-6)
})

