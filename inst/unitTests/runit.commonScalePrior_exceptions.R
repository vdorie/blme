cat("\n\nRUnit test cases for errors and warnings in blme:::parsePrior function, var.prior argument\n\n");


test.blme.parsePrior.cov.prior.exceptions <- function()
{
  generateData <- FALSE;
  
  testRoot <- file.path(.path.package(package="blme"), "unitTests");
  source(file.path(testRoot, "checkWarning.R"), TRUE);
  source(file.path(testRoot, "lmmData.R"), TRUE);
  
  options(warn = -1);
  lmerFit  <-  lmer(y ~ x.1 + (1 | g.1), control = list(maxIter = 0L));
  blmerFit <- blmer(y ~ x.1 + (1 | g.1), control = list(maxIter = 0L),
                    cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);
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
  
  checkException(parsePrior(blmerFit, var.prior = numeric(0)))
  checkException(parsePrior(blmerFit, var.prior = list(numeric(0))));
  checkException(parsePrior(blmerFit, var.prior = "not a prior"));
  
  checkException(parsePrior(blmerFit, var.prior = "point(()"));
  
  checkWarning(parsePrior(blmerFit, var.prior = "point(value = 2, notAParam = 0)"));
  checkException(parsePrior(blmerFit, var.prior = "point(value = 'not a number')"));
  checkException(parsePrior(blmerFit, var.prior = "point(value = 0)"));
  checkException(parsePrior(blmerFit, var.prior = "point(value = 2, posterior.scale = 'not a scale')"));

  checkException(parsePrior(blmerFit, var.prior = "inverse.gamma(-1)"));
  checkException(parsePrior(blmerFit, var.prior = "inverse.gamma(scale = -1)"));
  checkWarning(parsePrior(blmerFit, var.prior = "inverse.gamma(notAParam = 0)"));
  checkWarning(parsePrior(blmerFit, var.prior = "inverse.gamma(common.scale = 'anything')"));
  checkException(parsePrior(blmerFit, var.prior = "inverse.gamma(posterior.scale = 'not a scale')"));

  checkException(parsePrior(blmerFit, var.prior = "gamma(-1)"));
  checkException(parsePrior(blmerFit, var.prior = "gamma(rate = -1)"));
  checkWarning(parsePrior(blmerFit, var.prior = "gamma(notAParam = 0)"));
  checkException(parsePrior(blmerFit, var.prior = "gamma(posterior.scale = 'not a scale')"));
}
