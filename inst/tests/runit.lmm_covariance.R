cat("\n\nRUnit test cases for bmer:::blmer function with priors on the covariance of the modeled coefficients.\n\n");

test.blme.blmer.covarPrior <- function()
{
  source(system.file("common", "lmmData.R", package = "blme"));

  cov.prior <- "g.1 ~ gamma(rate = 0.5)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior);
  checkEquals(fit@pp$theta, 0.626025390625);

  cov.prior <- "g.1 ~ invgamma(scale = 2.0)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior);  
  checkEquals(fit@pp$theta, 0.939560546875);
  
  cov.prior <- "g.1 ~ wishart(scale = 2)";
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior);
  checkEquals(fit@pp$theta, c(0.677745102365688, -0.439777135132983, 1.48026251108622));

  cov.prior <- "g.1 ~ invwishart(scale = 2)";
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior);
  checkEquals(fit@pp$theta, c(0.627739008945695, -0.137563742254117, 1.05679359569432));

  cov.prior <- "g.1 ~ gamma(shape = 1.75, rate = 2, posterior.scale = 'var', common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior);
  checkEquals(fit@pp$theta, 0.435458984375);
  
  cov.prior <- "g.1 ~ invgamma(scale = 0.5, posterior.scale = 'var', common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior);
  checkEquals(fit@pp$theta, 0.460400390624999);

  cov.prior <- "g.1 ~ gamma(rate = 2, posterior.scale = 'sd', common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior);
  checkEquals(fit@pp$theta, 0.476779702210423);
  
  cov.prior <- "g.1 ~ invgamma(scale = 0.5, posterior.scale = 'sd', common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 | g.1), testData, cov.prior = cov.prior);
  checkEquals(fit@pp$theta, 0.452841796874999);

  cov.prior <- "g.1 ~ wishart(scale = 2, common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior);
  checkEquals(fit@pp$theta, c(0.63996739265564, -0.340538787006457, 1.34228986794088));
  
  cov.prior <- "g.1 ~ invwishart(scale = 2, common.scale = FALSE)";
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), testData, cov.prior = cov.prior);
  checkEquals(fit@pp$theta, c(0.505864621989816, -0.137623340382083, 0.979903012179649));
}
