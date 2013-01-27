cat("\n\nRUnit test cases for blme:::parsePrior function, var.prior argument\n\n");

test.blme.parsePrior.var.prior <- function()
{
  generateData <- FALSE;
  
  testRoot <- file.path(.path.package(package="blme"), "unitTests");
  source(file.path(testRoot, "lmmData.R"), TRUE);

  options(warn = -1);
  model1 <- blmer(y ~ x.1 + (1 | g.1), control = list(maxIter = 0L),
                  cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);

  parsePrior <- blme:::parsePrior;
  getScaleInt <- blme:::getScaleInt;
  getEnumOrder <- blme:::getEnumOrder;
  typeEnum <- blme:::typeEnum;
  familyEnum <- blme:::familyEnum;
  posteriorScaleEnum <- blme:::posteriorScaleEnum;
  commonScaleEnum <- blme:::commonScaleEnum;
  
  
  prior <- parsePrior(model1, var.prior = "point(2)");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 1);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "point"));
    checkEquals(prior@hyperparameters, 2);
  }

  prior <- parsePrior(model1, var.prior = "point(value = 3)");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 1);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "point"));
    checkEquals(prior@hyperparameters, 3);
  }

  prior <- parsePrior(model1, var.prior = "point(2, 'sd')");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 1);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "point"));
    checkEquals(prior@hyperparameters, 2);
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, "sd"),
                            getEnumOrder(commonScaleEnum, blme:::defaultCommonScaleCommonScale));
    checkEquals(prior@scales[1], scaleInt);
  }

  prior <- parsePrior(model1, var.prior = "point(2, posterior.scale = 'var')");
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 1);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "point"));
    checkEquals(prior@hyperparameters, 2);
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, "var"),
                            getEnumOrder(commonScaleEnum, blme:::defaultCommonScaleCommonScale));
    checkEquals(prior@scales[1], scaleInt);
  }

  prior <- parsePrior(model1, var.prior = "inverse.gamma")
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "inverse.gamma"));
  }

  prior <- parsePrior(model1, var.prior = "inverse.gamma(2, 0.5)")
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "inverse.gamma"));
    checkEquals(prior@hyperparameters, c(2, 0.5));
  }

  prior <- parsePrior(model1, var.prior = "inverse.gamma(posterior.scale = 'sd', shape = 3)")
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "inverse.gamma"));
    checkEquals(prior@hyperparameters, c(3, 0));
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, "sd"),
                            getEnumOrder(commonScaleEnum, blme:::defaultCommonScaleCommonScale));
    checkEquals(prior@scales[1], scaleInt);
  }

  prior <- parsePrior(model1, var.prior = "inverse.gamma")
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "inverse.gamma"));
  }

  prior <- parsePrior(model1, var.prior = "gamma(2, 0.5)")
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "gamma"));
    checkEquals(prior@hyperparameters, c(2, 0.5));
  }

  prior <- parsePrior(model1, var.prior = "gamma(posterior.scale = 'sd', shape = 1, rate = 0)")
  {
    checkEquals(prior@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior@families), 1);
    checkEquals(length(prior@scales), 1);
    checkEquals(length(prior@hyperparameters), 2);
    checkEquals(prior@families[1], getEnumOrder(familyEnum, "gamma"));
    checkEquals(prior@hyperparameters, c(1, 0));
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, "sd"),
                            getEnumOrder(commonScaleEnum, blme:::defaultCommonScaleCommonScale));
    checkEquals(prior@scales[1], scaleInt);
  }
}
