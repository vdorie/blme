cat("\n\nRUnit test cases for errors and warnings in blme:::parsePrior function, cov.prior argument\n\n");


test.blme.parsePrior.cov.prior.exceptions <- function()
{
  generateData <- FALSE;
  
  testRoot <- file.path(.path.package(package="blme"), "unitTests");
  source(file.path(testRoot, "checkWarning.R"), TRUE);
  source(file.path(testRoot, "lmmData.R"), TRUE);

  options(warn = -1);
  lmerFit  <-  lmer(y ~ x.1 + (1 | g.1), control = list(maxIter = 0L));
  blmerFit <- blmer(y ~ x.1 + (1 | g.1), control = list(maxIter = 0L),
                    cov.prior = NULL, fixef.prior = NULL);
  options(warn = 0);
  
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
  
  checkException(parsePrior(blmerFit, "correlation(()"));
  checkException(parsePrior(blmerFit, "g.1 ~ correlation(()"));
  checkException(parsePrior(blmerFit, "g.1 ~ correlation(1 ~ gamma()"));
  checkException(parsePrior(blmerFit, "g.1 ~ correlation(0 ~ gamma)"));
  checkException(parsePrior(blmerFit, "g.1 ~ correlation(notACoordinate ~ gamma)"));
  checkException(parsePrior(blmerFit, "g.1 ~ correlation(1 ~ wishart)"));
  
  checkWarning(parsePrior(blmerFit, "correlation(gamma(notAParam = 0))"))
  checkException(parsePrior(blmerFit, "correlation(gamma(shape = 'not a number'))"));
  checkException(parsePrior(blmerFit, "correlation(gamma(shape = -1))"));
  checkException(parsePrior(blmerFit, "correlation(gamma(rate = -1))"));
  checkException(parsePrior(blmerFit, "correlation(gamma(posterior.scale = 'kurtosis'))"));
  checkException(parsePrior(blmerFit, "correlation(gamma(common.scale = 'mu'))")); # mu in the Buddhist sense
  
  checkException(parsePrior(blmerFit, "spectral(()"));
  checkException(parsePrior(blmerFit, "g.1 ~ spectral(()"));
  checkException(parsePrior(blmerFit, "g.1 ~ spectral(1 ~ gamma()"));
  checkException(parsePrior(blmerFit, "g.1 ~ spectral(0 ~ gamma)"));
  
  checkWarning(parsePrior(blmerFit, "g.1 ~ spectral(gamma(notAParam = 0))"));
  
  checkException(parsePrior(blmerFit, "inverse.gamma(shape = 'not a number')"));
  checkException(parsePrior(blmerFit, "inverse.gamma(shape = -1)"));
  checkException(parsePrior(blmerFit, "inverse.gamma(scale = -1)"));
  
  checkException(parsePrior(blmerFit, "wishart(df = 'not a number')"));
  checkException(parsePrior(blmerFit, "wishart(df = 0)"));
  checkException(parsePrior(blmerFit, "wishart(scale = 'not a number')"));
  checkException(parsePrior(blmerFit, "wishart(scale = -0.01)"));
  checkException(parsePrior(blmerFit, "inverse.wishart(df = 'not a number')"));
  checkException(parsePrior(blmerFit, "inverse.wishart(df = 0)"));
  checkException(parsePrior(blmerFit, "inverse.wishart(inverse.scale = 'not a number')"));
  checkException(parsePrior(blmerFit, "inverse.wishart(inverse.scale = -0.01)"));

  checkException(parsePrior(blmerFit, "gamma(posterior.scale = 'not a scale')"));
  checkException(parsePrior(blmerFit, "gamma(common.scale = 'not a boolean')"));
  checkException(parsePrior(blmerFit, "gamma(data.scale = 'not a scale')"));


  options(warn = -1);
  blmerFit <- blmer(y ~ x.1 + (1 | g.1) + (0 + x.1 | g.1), control=list(maxIter = 0L));
  options(warn = 0);

  checkWarning(parsePrior(blmerFit, "g.1 ~ gamma"));
  checkException(parsePrior(blmerFit, "-1 ~ gamma"));
}
