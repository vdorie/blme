cat("\n\nRUnit test cases for bmer:::blmer function with priors on the covariance of the modeled coefficients.\n\n");

test.blme.blmer.covarPrior <- function()
{
  source(system.file("common", "lmmData.R", package = "blme"));

  control <- lmerControl(optimizer = "Nelder_Mead")
  lme4IsOld <- compareVersion(toString(packageVersion("lme4")), "1.1-4") < 0
  
  cov.prior <- "g.1 ~ gamma(rate = 0.5)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control);
  checkEquals(fit@pp$theta, 0.626025390625);

  cov.prior <- "g.1 ~ invgamma(scale = 2.0)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control);
  if (lme4IsOld) checkEquals(fit@theta, 0.939560546875)
  else checkEquals(fit@theta, 0.93955078125);
  
  cov.prior <- "g.1 ~ wishart(scale = 2)";
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control);
  checkEquals(fit@theta, c(0.677745102365688, -0.439777135132983, 1.48026251108622));

  cov.prior <- "g.1 ~ invwishart(scale = 2)";
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control);
  checkEquals(fit@theta, c(0.627739008945695, -0.137563742254117, 1.05679359569432));

  cov.prior <- "g.1 ~ gamma(shape = 1.75, rate = 2, posterior.scale = 'var', common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control);
  checkEquals(fit@theta, 0.435458984375);
  
  cov.prior <- "g.1 ~ invgamma(scale = 0.5, posterior.scale = 'var', common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control);
  if (lme4IsOld) checkEquals(fit@theta, 0.460400390624999)
  else checkEquals(fit@theta, 0.460410156249999);

  cov.prior <- "g.1 ~ gamma(rate = 2, posterior.scale = 'sd', common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control);
  checkEquals(fit@theta, 0.476779702210423);
  
  cov.prior <- "g.1 ~ invgamma(scale = 0.5, posterior.scale = 'sd', common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior, control = control);
  if (lme4IsOld) checkEquals(fit@theta, 0.452841796874999)
  else checkEquals(fit@theta, 0.452832031249999);

  cov.prior <- "g.1 ~ wishart(scale = 2, common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control);
  checkEquals(fit@pp$theta, c(0.63996739265564, -0.340538787006457, 1.34228986794088));
  
  cov.prior <- "g.1 ~ invwishart(scale = 2, common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior, control = control);
  checkEquals(fit@pp$theta, c(0.505864621989816, -0.137623340382083, 0.979903012179649));



  dwish <- function(R) {
    d <- nrow(R)
    nu <- d + 1 + 1.5;
    R.scale.inv <- diag(1e-2, d);
    
    const <- nu * (d * log(2) - 2 * sum(log(diag(R.scale.inv)))) +
      0.5 * d * (d - 1) * log(pi);
    for (i in 1:d) const <- const + 2 * lgamma(0.5 * (nu + 1.0 - i));
    
    det <- 2 * sum(log(diag(R)));
    
    const - (nu - d - 1) * det + sum((R %*% R.scale.inv)^2)
  }
  fit.prof <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, control = control,
                    cov.prior = wishart(scale = diag(1e4, q.k)));
  fit.cust <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, control = control,
                    cov.prior = custom(dwish, chol = TRUE, scale = "dev"));
  checkEquals(fit.prof@theta, fit.cust@theta, tolerance = 1e-6)
  
  fit.prof <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, control = control,
                    cov.prior = wishart(scale = diag(1e4, q.k), common.scale = FALSE));
  fit.cust <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, control = control,
                    cov.prior = custom(dwish, chol = TRUE, scale = "dev", common.scale = FALSE));
  checkEquals(c(fit.prof@pp$theta, fit.prof@devcomp$cmp[["sigmaREML"]]),
              c(fit.cust@pp$theta, fit.cust@devcomp$cmp[["sigmaREML"]]), tolerance = 5e-5)
}
