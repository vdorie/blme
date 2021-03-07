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
  } else {
    expect_equal(fit@theta, 1.26501385482573)
    expect_equal(fit@beta, c(0.873131216784199, 2.46647899095567, 1.32070496549675))
  }
})

test_that("bglmer fits test data with cov prior, matching previous version", {
  fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control,
                cov.prior = wishart)
  if (lme4Version < "1.1-4") {
    expect_equal(fit@theta, 2.96766525351892)
    expect_equal(fit@beta, c(1.0963789854971, 3.67790570859986, 1.75655010020603))
  } else {
    expect_equal(fit@theta, 2.96767284827046)
    expect_equal(fit@beta, c(1.0963789854971, 3.67790570859986, 1.75655010020603))
  }
})

test_that("bglmer runs with fixef horseshoe prior", {
  suppressWarnings(fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control,
                                 cov.prior = NULL, fixef.prior = horseshoe))
  
  expect_equal(fit@theta, 1.44188949480559)
  expect_equal(fit@beta, c(-3.00402517816968e-08, 2.38133218752372, 1.40658715296947))
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

