cat("\n\nRUnit test cases for bmer::blmer function with priors on the common scale/residual variance\n\n");

test.bmer.blmer.varPrior <- function()
{ 
  # eight schools
  y <- c(28, 8, -3, 7, -1, 1, 18, 12);
  sigma <- c(15, 10, 16, 11, 9, 11, 10, 18);

  y.z <- (y - mean(y)) / sigma;

  g <- 1:8;
  
  eightSchools <- blmer(y.z ~ 1 + (1 | g), resid.prior = point,
                        cov.prior = NULL, fixef.prior = NULL);

  checkEquals(eightSchools@pp$theta, 0);

  
  source(system.file("common", "lmmData.R", package = "blme"));
  
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2), testData,
               cov.prior = NULL, resid.prior = invgamma(2, 1.0));
  checkEquals(fit@pp$theta, c(0.720574821837846, -0.327950819407575, 1.48817902720885, 0.891726847920311, 0.314716339964349, -0.17288554050559, 1.02649696723723, 0.118403717814748, 0.205106180826464));
}
