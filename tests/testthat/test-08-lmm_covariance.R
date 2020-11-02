context("blmer numerical results with cov prior")

source(system.file("common", "lmmData.R", package = "blme"))
lme4Version <- packageVersion("lme4")
control <- lmerControl(optimizer = "Nelder_Mead")

test_that("blmer fits test data with gamma prior(TRUE), matching previous version", {
  cov.prior <- "g.1 ~ gamma(rate = 0.5)"
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control)
  if (lme4Version < "1.1-8") {
    expect_equal(fit@theta, 0.626025390625)
  } else {
    expect_equal(fit@theta, 0.626021723159)
  }
})

test_that("blmer fits test data with invgamma prior(TRUE), matching previous version", {
  cov.prior <- "g.1 ~ invgamma(scale = 2.0)"
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control)
  if (lme4Version < "1.1-4") {
    expect_equal(fit@theta, 0.93956054688)
  } else if (lme4Version < "1.1-8") {
    expect_equal(fit@theta, 0.93955078125)
  } else {
    expect_equal(fit@theta, 0.93955687941)
  }
})

test_that("blmer fits test data with wishart prior(TRUE), matching previous version", {
  cov.prior <- "g.1 ~ wishart(scale = 2)"
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control)
  expect_equal(fit@theta, c(0.677745102365688, -0.439777135132983, 1.48026251108622))
})

test_that("blmer fits test data with invwishart prior(TRUE), matching previous version", {
  cov.prior <- "g.1 ~ invwishart(scale = 2)"
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control)
  expect_equal(fit@theta, c(0.627739008945695, -0.137563742254117, 1.05679359569432))
})

test_that("blmer fits test data with gamma prior('var', FALSE),  matching previous version", {
  cov.prior <- "g.1 ~ gamma(shape = 1.75, rate = 2, posterior.scale = 'var', common.scale = FALSE)"
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control)
  if (lme4Version < "1.1-8") {
    expect_equal(fit@theta, 0.435458984375)
  } else {
    expect_equal(fit@theta, 0.435465082534)
  }
})

test_that("blmer fits test data with invgamma prior('var', FALSE), matching previous version", {
  cov.prior <- "g.1 ~ invgamma(scale = 0.5, posterior.scale = 'var', common.scale = FALSE)"
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control)
  if (lme4Version < "1.1-4") {
    expect_equal(fit@theta, 0.460400390625)
  } else if (lme4Version < "1.1-8") {
    expect_equal(fit@theta, 0.460410156250)
  } else {
    expect_equal(fit@theta, 0.460406488784)
  }
})

test_that("blmer fits test data with gamma prior('sd', FALSE), matching previous version", {
  cov.prior <- "g.1 ~ gamma(rate = 2, posterior.scale = 'sd', common.scale = FALSE)"
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control)
  if (lme4Version < "1.1-8") {
    expect_equal(fit@theta, 0.476779702210)
  } else {
    expect_equal(fit@theta, 0.476774614593)
  }
})

test_that("blmer fits test data with invgamma prior('sd', FALSE), matching previous version", {
  cov.prior <- "g.1 ~ invgamma(scale = 0.5, posterior.scale = 'sd', common.scale = FALSE)"
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control)
  if (lme4Version < "1.1-4") {
    expect_equal(fit@theta, 0.452841796875)
  } else if (lme4Version < "1.1-8") {
    expect_equal(fit@theta, 0.452832031250)
  } else {
    expect_equal(fit@theta, 0.452838129409)
  }
})

test_that("blmer fits test data with wishart prior(FALSE), matching previous version", {
  cov.prior <- "g.1 ~ wishart(scale = 2, common.scale = FALSE)"
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control)
  expect_equal(fit@pp$theta, c(0.63996739265564, -0.340538787006457, 1.34228986794088))
})

test_that("blmer fits test data with invwishart prior(FALSE), matching previous version", {
  cov.prior <- "g.1 ~ invwishart(scale = 2, common.scale = FALSE)"
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control)
  expect_equal(fit@pp$theta, c(0.505864621989816, -0.137623340382083, 0.979903012179649))
})

test_that("blmer fits test data with custom prior, matching builtin wishart", {
  dwish <- function(R) {
    d <- nrow(R)
    nu <- d + 1 + 1.5
    R.scale.inv <- diag(1e-2, d)
    
    const <- nu * (d * log(2) - 2 * sum(log(diag(R.scale.inv)))) +
      0.5 * d * (d - 1) * log(pi)
    for (i in 1:d) const <- const + 2 * lgamma(0.5 * (nu + 1.0 - i))
    
    det <- 2 * sum(log(diag(R)))
    
    const - (nu - d - 1) * det + sum((R %*% R.scale.inv)^2)
  }
  fit.prof <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, control = control,
                    cov.prior = wishart(scale = diag(1e4, q.k)))
  fit.cust <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, control = control,
                    cov.prior = custom(dwish, chol = TRUE, scale = "dev"))
  expect_equal(fit.prof@theta, fit.cust@theta, tolerance = 1e-6)
  
  fit.prof <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, control = control,
                    cov.prior = wishart(scale = diag(1e4, q.k), common.scale = FALSE))
  fit.cust <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, control = control,
                    cov.prior = custom(dwish, chol = TRUE, scale = "dev", common.scale = FALSE))
  expect_equal(c(fit.prof@pp$theta, fit.prof@devcomp$cmp[["sigmaREML"]]),
               c(fit.cust@pp$theta, fit.cust@devcomp$cmp[["sigmaREML"]]), tolerance = 5e-5)
})

test_that("blmer not confused by cov priors supplied in different forms", {
  control$optCtrl <- list(maxfun = 1L)
  control$checkConv <- NULL
  cov.prior <- list("g.1 ~ invgamma(scale = 2)", "g.2 ~ invgamma(scale = 3)")
  expect_is(suppressWarnings(blmer(y ~ x.1 + x.2 + (1 | g.1) + (1 | g.2), testData, cov.prior = cov.prior, control = control)),
            "blmerMod")
  cov.prior <- list(g.1 ~ invgamma(scale = 2), g.2 ~ invgamma(scale = 3))
  expect_is(suppressWarnings(blmer(y ~ x.1 + x.2 + (1 | g.1) + (1 | g.2), testData, cov.prior = cov.prior, control = control)),
            "blmerMod")
  expect_is(suppressWarnings(blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = wishart, control = control)), "blmerMod")
  expect_is(suppressWarnings(blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = wishart(), control = control)), "blmerMod")
})
