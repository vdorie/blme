cat("\n\nRUnit test cases for errors and warnings in blme:::parsePrior function, fixef.prior argument\n\n");


test.blme.parsePrior.cov.prior.exceptions <- function()
{
  source(system.file("common", "lmmData.R", package = "blme"));

  fit <- blmer(y ~ x.1 + (1 | g.1), testData,
               cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL);
  
  RUnitOptions <- getOption("RUnit");
  RUnitOptions$silent <- TRUE;
  options("RUnit" = RUnitOptions);

  parsePrior <- blme:::parsePrior;
  
  checkException(parsePrior(fit, fixef.prior = "normal(common.scale = 'crazy')"));
  checkException(parsePrior(fit, fixef.prior = "normal(cov = diag(3))"));
  negDefiniteMatrix <- matrix(c(1, 0, 0, -0.1), 2, 2);
  checkException(parsePrior(fit, fixef.prior = "normal(cov = negDefiniteMatrix)"));
  asymmetricMatrix <- matrix(c(1, 0.5, 0.3, 0.7), 2, 2);
  checkException(parsePrior(fit, fixef.prior = "normal(cov = asymmetricMatrix)"));

  checkException(parsePrior(fit, fixef.prior = "t"))

  fit <- blmer(y ~ x.1 + (1 | g.1), testData, REML = FALSE,
               cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL);
  checkException(parsePrior(fit, fixef.prior = "t(df = 0)"));
  checkException(parsePrior(fit, fixef.prior = "t(scale = c(-1, 2))"));
}
