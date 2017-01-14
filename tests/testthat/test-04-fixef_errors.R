context("b(g)lmer, fixef.prior argument")

test_that("fixef.prior argument raises appropriate errors for blmer fits", {
  source(system.file("common", "lmmData.R", package = "blme"))

  fit <- blmer(y ~ x.1 + (1 | g.1), testData,
               cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL)
  
  parsePrior <- blme:::parsePrior
  
  expect_error(parsePrior(fit, fixef.prior = "normal(common.scale = 'crazy')"))
  expect_error(parsePrior(fit, fixef.prior = "normal(cov = diag(3))"))
  negDefiniteMatrix <- matrix(c(1, 0, 0, -0.1), 2, 2)
  expect_error(parsePrior(fit, fixef.prior = "normal(cov = negDefiniteMatrix)"))
  asymmetricMatrix <- matrix(c(1, 0.5, 0.3, 0.7), 2, 2)
  expect_error(parsePrior(fit, fixef.prior = "normal(cov = asymmetricMatrix)"))

  expect_error(parsePrior(fit, fixef.prior = "t"))

  fit <- blmer(y ~ x.1 + (1 | g.1), testData, REML = FALSE,
               cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL)
  expect_error(parsePrior(fit, fixef.prior = "t(df = 0)"))
  expect_error(parsePrior(fit, fixef.prior = "t(scale = c(-1, 2))"))
  expect_error(parsePrior(fit, fixef.prior = "horseshoe(global.shrinkage = -1)"))
  
  expect_error(blmer(y ~ x.1 + (1 | g.1), testData, REML = TRUE,
               cov.prior = NULL, fixef.prior = t, resid.prior = NULL))
  expect_error(blmer(y ~ x.1 + (1 | g.1), testData, REML = TRUE,
               cov.prior = NULL, fixef.prior = horseshoe(), resid.prior = NULL))
})

test_that("fixef.prior argument raises appropriate errors for bglmer fits", {
  source(system.file("common", "glmmData.R", package = "blme"))
  
  parsePrior <- blme:::parsePrior
  
  bglmerFit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(),
                      cov.prior = NULL)
  
  expect_error(parsePrior(bglmerFit, fixef.prior = normal(common.scale = TRUE)))
  expect_error(parsePrior(bglmerFit, fixef.prior = t(common.scale = TRUE)))
  expect_error(parsePrior(bglmerFit, fixef.prior = horseshoe(common.scale = TRUE)))
})
