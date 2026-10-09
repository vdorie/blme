## for old versions of lme4, the refit for g/lmer doesn't move very
## far from the fit so we just suppress the test

lme4Version <- packageVersion("lme4")
if (lme4Version >= "1.1-6") {
  context("refit generic for blmerMod and bglmerMod classes")

  source(system.file("common", "lmmData.R", package = "blme"))
  control <- lmerControl(optimizer = "Nelder_Mead")
  
  cov.prior <- "g.1 ~ wishart(scale = 2)"
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control)
  
  test_that("refit for blmer matches original, not lmer", {
    blmerRefit <- refit(fit)
    lmerRefit  <- getS3method("refit", "merMod")(fit)
    
    expect_equal(fit@theta, blmerRefit@theta, tolerance = 1.0e-02)
    expect_equal(fit@beta,  blmerRefit@beta,  tolerance = 1.0e-03)
    
    expect_false(all(abs(fit@theta - lmerRefit@theta) <= 1.0e-02))
    expect_false(all(abs(fit@beta  - lmerRefit@beta)  <= 1.0e-03))
    
    ## REML adjustment depends on the number of fixed effects
    expect_equal(blmerRefit@devcomp$cmp[["REML"]], fit@devcomp$cmp[["REML"]], tolerance = 1.0e-4)
  })

  test_that("ml refit matches ml fit", {
    mlRefit <- refitML(fit)
    mlFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control, REML = FALSE)

    expect_equal(mlRefit@theta, mlFit@theta, tolerance = 1.0e-02)
    expect_equal(mlRefit@beta,  mlFit@beta,  tolerance = 1.0e-03)

    expect_false(all(abs(fit@theta - mlRefit@theta) <= 1.0e-02))
    expect_false(all(abs(fit@beta  - mlFit@beta)  <= 1.0e-03))
  })

  test_that("ml refit takes refit arguments and numerically optimized sigma", {
    fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control,
                 fixef.prior = normal(common.scale = FALSE))
    mlRefit <- refitML(fit, maxit = 50L)
    mlFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control,
                   fixef.prior = normal(common.scale = FALSE), REML = FALSE)
    
    expect_false(is.na(logLik(mlRefit)))
    expect_equal(mlRefit@theta, mlFit@theta, tolerance = 1.0e-02)
    expect_equal(sigma(mlRefit), sigma(mlFit), tolerance = 1.0e-03)
    
    ## control goes to the optimizer, as in lme4
    mlRefit <- suppressWarnings(refitML(fit, control = list(maxfun = 5L)))
    expect_equal(mlRefit@optinfo$feval, 5L)
  })
  
  test_that("refit for bglmer with nAGQ = 0 matches original", {
    source(system.file("common", "glmmData.R", package = "blme"))
    
    fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), nAGQ = 0)
    bglmerRefit <- refit(fit)
    
    expect_equal(fit@theta, bglmerRefit@theta, tolerance = 1.0e-2)
    expect_equal(fit@beta,  bglmerRefit@beta,  tolerance = 1.0e-3)
  })
  
  test_that("refit for bglmer matches original, not glmer", {
    source(system.file("common", "glmmData.R", package = "blme"))
    control <- if (lme4Version >= "1.1-8")
                 glmerControl(optimizer = "Nelder_Mead", nAGQ0initStep = FALSE)
               else
                 glmerControl(optimizer = "Nelder_Mead")
    
    fit <- bglmer(y ~ x.1 + x.2 + (1 | g), testData, family = binomial(), control = control,
                  cov.prior = wishart)
    
    bglmerRefit <- refit(fit)
    glmerRefit  <- getS3method("refit", "merMod")(fit)
    
    expect_equal(fit@theta, bglmerRefit@theta, tolerance = 1.0e-3)
    expect_equal(fit@beta,  bglmerRefit@beta,  tolerance = 1.0e-3)
    
    expect_false(all(abs(fit@theta - glmerRefit@theta) <= 1.0e-3))
    expect_false(all(abs(fit@beta  - glmerRefit@beta)  <= 1.0e-3))
  })
}
