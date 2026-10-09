context("glmer and bglmer")

test_that("bglmer matches glmer exactly", {
  source(system.file("common", "glmmData.R", package = "blme"))
  
  control <- glmerControl(optimizer = "Nelder_Mead")
  
   glmerFit <-  glmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control)
  bglmerFit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control,
                      cov.prior = NULL)
  
  expect_equal(glmerFit@theta, bglmerFit@theta)
  expect_equal(glmerFit@beta, bglmerFit@beta)
  expect_equal(glmerFit@u, bglmerFit@u)
})

test_that("bglmer fits with nAGQ = 0 and with starting values when skipping that step", {
  source(system.file("common", "glmmData.R", package = "blme"))
  
   glmerFit <-  glmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), nAGQ = 0)
  bglmerFit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), nAGQ = 0,
                      cov.prior = NULL)
  expect_equal(glmerFit@theta, bglmerFit@theta)
  expect_equal(glmerFit@beta, bglmerFit@beta)
  fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), nAGQ = 0)
  expect_is(fit, "bglmerMod")
  expect_length(fit@optinfo$val, length(fit@theta))
  expect_error(bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), nAGQ = 0,
                      fixef.prior = normal), "nAGQ")
  
  if (packageVersion("lme4") >= "1.1-8") {
    fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(),
                  control = glmerControl(nAGQ0initStep = FALSE),
                  start = list(theta = 1, fixef = c(0, 0, 0)))
    expect_is(fit, "bglmerMod")
  }
})

test_that("bglmer doesn't report a singular fit away from the boundary", {
  source(system.file("common", "glmmData.R", package = "blme"))
  
  expect_message(fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial()), NA)
  expect_false(isSingular(fit))
})
