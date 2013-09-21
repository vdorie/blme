cat("\n\nRUnit test cases for errors and warnings in blme:::parsePrior function, var.prior argument\n\n");


test.blme.parsePrior.cov.prior.exceptions <- function()
{
  source(system.file("common", "lmmData.R", package = "blme"));
  source(system.file("common", "checkWarning.R", package = "blme"));

  lmerFit  <-  lmer(y ~ x.1 + (1 | g.1), testData);
  blmerFit <- blmer(y ~ x.1 + (1 | g.1), testData,
                    cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL);
  
  RUnitOptions <- getOption("RUnit");
  RUnitOptions$silent <- TRUE;
  options("RUnit" = RUnitOptions);

  parsePrior <- blme:::parsePrior;
  
  checkException(parsePrior(blmerFit, resid.prior = numeric(0)));
  checkException(parsePrior(blmerFit, resid.prior = list(numeric(0))));
  checkException(parsePrior(blmerFit, resid.prior = "not a prior"));
  
  checkException(parsePrior(blmerFit, resid.prior = "point(()"));
  
  checkException(parsePrior(blmerFit, resid.prior = "point(value = 2, notAParam = 0)"));
  checkException(parsePrior(blmerFit, resid.prior = "point(value = 'not a number')"));
  checkException(parsePrior(blmerFit, resid.prior = "point(value = 0)"));
  checkException(parsePrior(blmerFit, resid.prior = "point(value = 2, posterior.scale = 'not a scale')"));

  checkException(parsePrior(blmerFit, resid.prior = "invgamma(-1)"));
  checkException(parsePrior(blmerFit, resid.prior = "invgamma(scale = -1)"));
  checkException(parsePrior(blmerFit, resid.prior = "invgamma(notAParam = 0)"));
  checkException(parsePrior(blmerFit, resid.prior = "invgamma(common.scale = 'anything')"));
  checkException(parsePrior(blmerFit, resid.prior = "invgamma(posterior.scale = 'not a scale')"));

  checkException(parsePrior(blmerFit, resid.prior = "gamma(-1)"));
  checkException(parsePrior(blmerFit, resid.prior = "gamma(rate = -1)"));
  checkException(parsePrior(blmerFit, resid.prior = "gamma(notAParam = 0)"));
  checkException(parsePrior(blmerFit, resid.prior = "gamma(posterior.scale = 'not a scale')"));
}
