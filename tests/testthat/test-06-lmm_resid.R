context("blmer numerical results with residual variance prior")

test_that("blmer fits the eight schools example correctly", {
  # eight schools
  y <- c(28, 8, -3, 7, -1, 1, 18, 12)
  sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)

  y.z <- (y - mean(y)) / sigma

  g <- 1:8
  
  eightSchools <- blmer(y.z ~ 1 + (1 | g), resid.prior = point,
                        cov.prior = NULL, fixef.prior = NULL)

  expect_equal(eightSchools@theta, 0, tolerance = 1.0e-7)
})

test_that("blmer fits test data with inv.gamma prior, matching previous version", {
  source(system.file("common", "lmmData.R", package = "blme"))
  
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2), testData, control = lmerControl(optimizer = "bobyqa"),
               cov.prior = NULL, resid.prior = invgamma(2, 1.0))
  expect_equal(fit@pp$theta, c(0.725321928185923, -0.251272308917427, 1.5828609906233, 0.946932542474828, 0.467970716580088, -0.183212783510381, 1.07158442297183, 0.122067368879505, 0.223238050522642),
               tolerance = 1.0e-5)
})
