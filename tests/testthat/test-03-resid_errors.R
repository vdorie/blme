context("blmer, resid.prior argument")

test_that("resid.prior argument raises apprprioate errors", {
  source(system.file("common", "lmmData.R", package = "blme"))

  blmerFit <- blmer(y ~ x.1 + (1 | g.1), testData,
                    cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL)
  
  parsePrior <- blme:::parsePrior
  
  expect_error(parsePrior(blmerFit, resid.prior = numeric(0)))
  expect_error(parsePrior(blmerFit, resid.prior = list(numeric(0))))
  expect_error(parsePrior(blmerFit, resid.prior = "not a prior"))
  
  expect_error(parsePrior(blmerFit, resid.prior = "point(()"))
  
  expect_error(parsePrior(blmerFit, resid.prior = "point(value = 2, notAParam = 0)"))
  expect_error(parsePrior(blmerFit, resid.prior = "point(value = 'not a number')"))
  expect_error(parsePrior(blmerFit, resid.prior = "point(value = 0)"))
  expect_error(parsePrior(blmerFit, resid.prior = "point(value = 2, posterior.scale = 'not a scale')"))

  expect_error(parsePrior(blmerFit, resid.prior = "invgamma(-1)"))
  expect_error(parsePrior(blmerFit, resid.prior = "invgamma(scale = -1)"))
  expect_error(parsePrior(blmerFit, resid.prior = "invgamma(notAParam = 0)"))
  expect_error(parsePrior(blmerFit, resid.prior = "invgamma(common.scale = 'anything')"))
  expect_error(parsePrior(blmerFit, resid.prior = "invgamma(posterior.scale = 'not a scale')"))

  expect_error(parsePrior(blmerFit, resid.prior = "gamma(-1)"))
  expect_error(parsePrior(blmerFit, resid.prior = "gamma(rate = -1)"))
  expect_error(parsePrior(blmerFit, resid.prior = "gamma(notAParam = 0)"))
  expect_error(parsePrior(blmerFit, resid.prior = "gamma(posterior.scale = 'not a scale')"))
})
