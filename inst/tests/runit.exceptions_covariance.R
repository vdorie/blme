cat("\n\nRUnit test cases for errors and warnings in blme:::parsePrior function, cov.prior argument\n\n");


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
  
  # Morally speaking, parsePrior isn't exposed to the user
  # so perhaps this first set of tests is excessive.
  checkException(parsePrior());
  checkException(parsePrior(NULL));
  checkException(parsePrior(notAValidObject));
  checkException(parsePrior(lmerFit));
  
  checkException(parsePrior(blmerFit, numeric(0)))
  checkException(parsePrior(blmerFit, list(numeric(0))));
  checkException(parsePrior(blmerFit, "not a prior"));
  checkException(parsePrior(blmerFit, list("not a", "prior")));
  checkException(parsePrior(blmerFit, "notAGroup ~ gamma"));
  
  checkException(parsePrior(blmerFit, "invgamma(shape = 'not a number')"));
  checkException(parsePrior(blmerFit, "invgamma(shape = -1)"));
  checkException(parsePrior(blmerFit, "invgamma(scale = -1)"));
  
  checkException(parsePrior(blmerFit, "wishart(df = 'not a number')"));
  checkException(parsePrior(blmerFit, "wishart(df = 0)"));
  checkException(parsePrior(blmerFit, "wishart(scale = 'not a number')"));
  checkException(parsePrior(blmerFit, "wishart(scale = -0.01)"));
  checkException(parsePrior(blmerFit, "invwishart(df = 'not a number')"));
  checkException(parsePrior(blmerFit, "invwishart(df = 0)"));
  checkException(parsePrior(blmerFit, "invwishart(scale = 'not a number')"));
  checkException(parsePrior(blmerFit, "invwishart(scale = -0.01)"));

  checkException(parsePrior(blmerFit, "gamma(posterior.scale = 'not a scale')"));
  checkException(parsePrior(blmerFit, "gamma(common.scale = 'not a boolean')"));

  blmerFit <- blmer(y ~ x.1 + (1 + x.1 | g.1), testData,
                    cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL);
  checkWarning(parsePrior(blmerFit, "gamma"));
}
