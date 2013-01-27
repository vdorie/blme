cat("\n\nRUnit test cases for blme:::parsePrior function, fixef.prior argument\n\n");


test.blme.parsePrior.fixef.prior <- function()
{
  generateData <- FALSE;
  
  testRoot <- file.path(.path.package(package="blme"), "unitTests");
  
  source(file.path(testRoot, "lmmData.R"), TRUE);
  source(file.path(testRoot, "checkWarning.R"), TRUE);

  # ignore the false convergence
  options(warn = -1);
  model1 <- blmer(y ~ x.1 + (1 | g.1), control = list(maxIter = 0L),
                  cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);

  parsePrior <- blme:::parsePrior;
  getEnumOrder <- blme:::getEnumOrder;
  getScaleInt <- blme:::getScaleInt;
  typeEnum <- blme:::typeEnum;
  familyEnum <- blme:::familyEnum;
  posteriorScaleEnum <- blme:::posteriorScaleEnum;
  commonScaleEnum <- blme:::commonScaleEnum;
  
  prior <- parsePrior(model1, fixef.prior = "normal(data.scale = 'absolute')");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 3);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "normal"));
  }

  prior <- parsePrior(model1, fixef.prior = "normal(1.5)");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "normal"));
    checkEquals(prior@hyperparameters, c(2.0 * 2.0 * log(1.5), 1 / 1.5));
  }

  prior <- parsePrior(model1, fixef.prior = "normal(1.5, 'false')");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "normal"));
    checkEquals(prior@hyperparameters, c(2.0 * 2.0 * log(1.5), 1 / 1.5));
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, blme:::defaultUnmodeledCoefficientPosteriorScale),
                            getEnumOrder(commonScaleEnum, 'false'));
    checkEquals(prior@scales[1], scaleInt);
  }

  prior <- parsePrior(model1, fixef.prior = "normal(diag(2), common.scale = 'false')");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 1 + 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "normal"));
    checkEquals(prior@hyperparameters, c(0, rep(1, 2)));
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, blme:::defaultUnmodeledCoefficientPosteriorScale),
                            getEnumOrder(commonScaleEnum, 'false'));
    checkEquals(prior@scales[1], scaleInt);
  }

  testCovariance <- matrix(c(2.51370232565668, 1.41538355022191, 1.41538355022191, 2.91044687425429), 2, 2);
  prior <- parsePrior(model1, fixef.prior = "normal(cov = testCovariance, common.scale = 'false')");
  {
    leftFactor <- solve(t(chol(testCovariance)));
    covariance.inv <- solve(testCovariance);
    logDetCov <- as.numeric(determinant(testCovariance, logarithm=TRUE)$modulus);
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 1 + 2 * 4);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "normal"));
    checkEquals(prior@hyperparameters, c(logDetCov, leftFactor, covariance.inv));
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, blme:::defaultUnmodeledCoefficientPosteriorScale),
                            getEnumOrder(commonScaleEnum, 'false'));
    checkEquals(prior@scales[1], scaleInt);
  }
  
  prior <- parsePrior(model1, fixef.prior = "normal(cov = c(1.1, 0.9))");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 1 + 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "normal"));
    checkEquals(prior@hyperparameters, c(log(1.1) + log(0.9), 1 / sqrt(c(1.1, 0.9))));
  }

  prior <- parsePrior(model1, fixef.prior = "normal(sd = 1.2, data.scale = 'free')");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 1 + 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "normal"));
    checkEquals(prior@hyperparameters, c(2 * 2 * log(1.2) + 2 * log(var(y)) - log(var(x.1)),
                                         1 / c(1.2 * sd(y), 1.2 * sd(y) / sd(x.1))));
  }
  
  # exceptions and warnings
  RUnitOptions <- getOption("RUnit");
  RUnitOptions$silent <- TRUE;
  options("RUnit" = RUnitOptions);
  
  checkException(blme:::parsePrior(model1, fixef.prior = "normal(common.scale = 'crazy')"));
  checkException(blme:::parsePrior(model1, fixef.prior = "normal(cov = diag(3))"));
  negDefiniteMatrix <- matrix(c(1, 0, 0, -0.1), 2, 2);
  checkException(blme:::parsePrior(model1, fixef.prior = "normal(cov = negDefiniteMatrix)"));
  asymmetricMatrix <- matrix(c(1, 0.5, 0.3, 0.7), 2, 2);
  checkException(blme:::parsePrior(model1, fixef.prior = "normal(cov = asymmetricMatrix)"));
}
