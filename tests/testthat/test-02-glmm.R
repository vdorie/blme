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
