context("blmer, cov.prior argument")


test_that("cov.prior argument raises appropriate errors", {
  source(system.file("common", "lmmData.R", package = "blme"))

  lmerFit  <-  lmer(y ~ x.1 + (1 | g.1), testData)
  blmerFit <- blmer(y ~ x.1 + (1 | g.1), testData,
                    cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL)
  
  parsePrior <- blme:::parsePrior
  
  # Morally speaking, parsePrior isn't exposed to the user
  # so perhaps this first set of tests is excessive.
  expect_error(parsePrior())
  expect_error(parsePrior(NULL))
  expect_error(parsePrior(notAValidObject))
  expect_error(parsePrior(lmerFit))
  
  expect_error(parsePrior(blmerFit, numeric(0)))
  expect_error(parsePrior(blmerFit, list(numeric(0))))
  expect_error(parsePrior(blmerFit, "not a prior"))
  expect_error(parsePrior(blmerFit, list("not a", "prior")))
  expect_error(parsePrior(blmerFit, "notAGroup ~ gamma"))
  
  expect_error(parsePrior(blmerFit, "invgamma(shape = 'not a number')"))
  expect_error(parsePrior(blmerFit, "invgamma(shape = -1)"))
  expect_error(parsePrior(blmerFit, "invgamma(scale = -1)"))
  
  expect_error(parsePrior(blmerFit, "wishart(df = 'not a number')"))
  expect_error(parsePrior(blmerFit, "wishart(df = 0)"))
  expect_error(parsePrior(blmerFit, "wishart(scale = 'not a number')"))
  expect_error(parsePrior(blmerFit, "wishart(scale = -0.01)"))
  expect_error(parsePrior(blmerFit, "invwishart(df = 'not a number')"))
  expect_error(parsePrior(blmerFit, "invwishart(df = 0)"))
  expect_error(parsePrior(blmerFit, "invwishart(scale = 'not a number')"))
  expect_error(parsePrior(blmerFit, "invwishart(scale = -0.01)"))

  expect_error(parsePrior(blmerFit, "gamma(posterior.scale = 'not a scale')"))
  expect_error(parsePrior(blmerFit, "gamma(common.scale = 'not a boolean')"))

  blmerFit <- blmer(y ~ x.1 + (1 + x.1 | g.1), testData,
                    cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL)
  expect_warning(parsePrior(blmerFit, "gamma"))
})
