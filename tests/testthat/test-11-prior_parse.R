library(blme)
source(system.file("common", "lmmData.R", package = "blme"))
lme4Version <- packageVersion("lme4")
control <- lmerControl(optimizer = "Nelder_Mead")
control$optCtrl <- list(maxfun = 1L)
control$checkConv <- NULL

test_that("blmer finds prior options specified as variables in global env", {
  g1_rate <- 0.5
  cov.prior <- "g.1 ~ gamma(rate = g1_rate)"
  expect_is(
    suppressWarnings(blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control)),
    "blmerMod"
  )

  expect_is(
    suppressWarnings(blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = g.1 ~ gamma(rate = g1_rate), control = control)),
    "blmerMod"
  )
})

test_that("blmer finds prior options when fit in a function", {
  fit_fn <- function() {
    g1_rate <- 0.5
    cov.prior <- "g.1 ~ gamma(rate = g1_rate)"
    blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control)
  }
  expect_is(suppressWarnings(fit_fn()), "blmerMod")

  fit_fn <- function() {
    g1_rate <- 0.5
    blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = g.1 ~ gamma(rate = g1_rate), control = control)
  }
  expect_is(suppressWarnings(fit_fn()), "blmerMod")

  g1_rate <- 0.5
  fit_fn <- function() {
    blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = g.1 ~ gamma(rate = g1_rate), control = control)
  }
  expect_is(suppressWarnings(fit_fn()), "blmerMod")
})

# test thanks to Jacob Grytzka
test_that("blmer finds prior options when fit in a function, horseshoe prior specifically", {
  skip_if_not_installed("expint")
  fit_fn <- function() {
    pen <- 5
    blmer(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + (1|Species),
          data = iris, REML = FALSE,
          fixef.prior = horseshoe(mean = 1e-5,
                                  global.shrinkage = pen,
                                  common.scale = TRUE)
    )
  }
  expect_is(suppressWarnings(fit_fn()), "blmerMod")

  pen <- 5
  fit_fn <- function() {
    blmer(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + (1|Species),
          data = iris, REML = FALSE,
          fixef.prior = horseshoe(mean = 1e-5,
                                  global.shrinkage = pen,
                                  common.scale = TRUE)
    )
  }
  expect_is(suppressWarnings(fit_fn()), "blmerMod")
})

