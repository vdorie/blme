cat("\n\nRUnit test cases for blme:::parsePrior function, cov.prior argument\n\n");

test.blme.parsePrior.cov.prior <- function()
{
  generateData <- FALSE;
  
  testRoot <- file.path(.path.package(package="blme"), "unitTests");
  source(file.path(testRoot, "lmmData.R"), TRUE);

  # snag some non-exported stuff from blme to make cleaner
  parsePrior <- blme:::parsePrior;
  getEnumOrder <- blme:::getEnumOrder;
  getScaleInt <- blme:::getScaleInt;
  typeEnum <- blme:::typeEnum;
  familyEnum <- blme:::familyEnum;
  posteriorScaleEnum <- blme:::posteriorScaleEnum;
  commonScaleEnum <- blme:::commonScaleEnum;
  
  options(warn = -1);
  model1 <- blmer(y ~ x.1 + (1 | g.1), control = list(maxIter = 0L),
                  cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);
  
  prior <- parsePrior(model1, cov.prior = "gamma");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "gamma"));
  }

  prior <- parsePrior(model1, cov.prior = "gamma(3)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "gamma"));
    checkEquals(prior[[1]]@hyperparameters[1], 3);
  }
  
  prior <- parsePrior(model1, cov.prior = "gamma(2, 1.8)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "gamma"));
    checkEquals(prior[[1]]@hyperparameters, c(2, 1.8));
  }
  
  prior <- parsePrior(model1, cov.prior = "gamma(rate = 1.8, shape = 2, posterior.scale = 'var')")
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2);
    checkEquals(prior[[1]]@families[1], getEnumOrder(blme:::familyEnum, "gamma"));
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, "var"),
                            getEnumOrder(commonScaleEnum, blme:::defaultDirectCommonScale));
    checkEquals(prior[[1]]@scales[1], scaleInt);
    checkEquals(prior[[1]]@hyperparameters, c(2, 1.8));
  }

  prior <- parsePrior(model1, cov.prior = "gamma(common.scale = 'false')")
  {
    checkEquals(prior[[1]]@type, getEnumOrder(blme:::typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "gamma"));
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, blme:::defaultDirectPosteriorScale),
                            getEnumOrder(commonScaleEnum, "false"));
    checkEquals(prior[[1]]@scales, scaleInt);
  }

  prior <- parsePrior(model1, cov.prior = "gamma(rate = 0)")
  {
    checkEquals(prior[[1]]@type, getEnumOrder(blme:::typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "gamma"));
    checkEquals(prior[[1]]@hyperparameters[2], 0);
  }

  prior <- parsePrior(model1, cov.prior = "inverse.gamma(scale = 0)")
  {
    checkEquals(prior[[1]]@type, getEnumOrder(blme:::typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "inverse.gamma"));
    checkEquals(prior[[1]]@hyperparameters[2], 0);
  }
  
  options(warn = -1);
  model2 <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                  control = list(maxIter = 0L),
                  cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);
  d.1 <- 2;
  d.2 <- 3;
  
  prior <- parsePrior(model2, cov.prior = "g.1 ~ spectral, g.2 ~ none");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "spectral"));
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "none"));
    checkEquals(length(prior[[1]]@families), d.1);
    checkEquals(length(prior[[1]]@scales), d.1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 * d.1);
    checkEquals(length(prior[[2]]@families), 0);
    checkEquals(length(prior[[2]]@scales), 0);
    checkEquals(length(prior[[2]]@hyperparameters), 0);
  }

  prior <- parsePrior(model2, cov.prior = "g.1 ~ spectral(inverse.gamma), g.2 ~ none");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "spectral"));
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "none"));
    checkEquals(length(prior[[1]]@families), d.1);
    checkEquals(length(prior[[1]]@scales), d.1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 * d.1);
    checkEquals(length(prior[[2]]@families), 0);
    checkEquals(length(prior[[2]]@scales), 0);
    checkEquals(length(prior[[2]]@hyperparameters), 0);

    checkEquals(prior[[1]]@families, rep(getEnumOrder(familyEnum, "inverse.gamma"), d.1));
  }

  prior <- parsePrior(model2, cov.prior = "g.1 ~ spectral(gamma(shape = 3, rate = .25)), g.2 ~ none");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "spectral"));
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "none"));
    checkEquals(length(prior[[1]]@families), d.1);
    checkEquals(length(prior[[1]]@scales), d.1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 * d.1);
    checkEquals(length(prior[[2]]@families), 0);
    checkEquals(length(prior[[2]]@scales), 0);
    checkEquals(length(prior[[2]]@hyperparameters), 0);

    checkEquals(prior[[1]]@families, rep(getEnumOrder(familyEnum, "gamma"), d.1));
    checkEquals(prior[[1]]@hyperparameters, rep(c(3, 0.25), d.1));
  }

  prior <- parsePrior(model2, cov.prior = "g.1 ~ spectral(gamma(shape = 3, rate = .25, posterior.scale = 'sd')), g.2 ~ none");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "spectral"));
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "none"));
    checkEquals(length(prior[[1]]@families), d.1);
    checkEquals(length(prior[[1]]@scales), d.1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 * d.1);
    checkEquals(length(prior[[2]]@families), 0);
    checkEquals(length(prior[[2]]@scales), 0);
    checkEquals(length(prior[[2]]@hyperparameters), 0);

    checkEquals(prior[[1]]@families, rep(getEnumOrder(familyEnum, "gamma"), d.1));
    checkEquals(prior[[1]]@hyperparameters, rep(c(3, 0.25), d.1));
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, "sd"),
                            getEnumOrder(commonScaleEnum, blme:::defaultDirectCommonScale));
    checkEquals(prior[[1]]@scales, rep(scaleInt, d.1));
  }
  
  prior <- parsePrior(model2, cov.prior = "spectral(gamma(shape = 3, rate = .25, posterior.scale = 'sd'))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "spectral"));
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "spectral"));
    checkEquals(length(prior[[1]]@families), d.1);
    checkEquals(length(prior[[1]]@scales), d.1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 * d.1);
    checkEquals(length(prior[[2]]@families), d.2);
    checkEquals(length(prior[[2]]@scales), d.2);
    checkEquals(length(prior[[2]]@hyperparameters), 2 * d.2);

    
    checkEquals(prior[[1]]@families, rep(getEnumOrder(familyEnum, "gamma"), d.1));
    checkEquals(prior[[1]]@hyperparameters, rep(c(3, 0.25), d.1));

    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, "sd"),
                            getEnumOrder(commonScaleEnum, blme:::defaultDirectCommonScale));
    checkEquals(prior[[1]]@scales, rep(scaleInt, d.1));

    
    checkEquals(prior[[2]]@families, rep(getEnumOrder(familyEnum, "gamma"), d.2));
    checkEquals(prior[[2]]@hyperparameters, rep(c(3, 0.25), d.2));

    checkEquals(prior[[2]]@scales, rep(scaleInt, d.2));
  }

  prior <- parsePrior(model2, cov.prior = "spectral(gamma(shape = 3, rate = .25, common.scale = 'false'))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "spectral"));
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "spectral"));
    checkEquals(length(prior[[1]]@families), d.1);
    checkEquals(length(prior[[1]]@scales), d.1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 * d.1);
    checkEquals(length(prior[[2]]@families), d.2);
    checkEquals(length(prior[[2]]@scales), d.2);
    checkEquals(length(prior[[2]]@hyperparameters), 2 * d.2);

    
    checkEquals(prior[[1]]@families, rep(getEnumOrder(familyEnum, "gamma"), d.1));
    checkEquals(prior[[1]]@hyperparameters, rep(c(3, 0.25), d.1));

    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, blme:::defaultDirectPosteriorScale),
                            getEnumOrder(commonScaleEnum, 'false'));
    checkEquals(prior[[1]]@scales, rep(scaleInt, d.1));

    
    checkEquals(prior[[2]]@families, rep(getEnumOrder(familyEnum, "gamma"), d.2));
    checkEquals(prior[[2]]@hyperparameters, rep(c(3, 0.25), d.2));

    checkEquals(prior[[2]]@scales, rep(scaleInt, d.2));
  }

  prior <- parsePrior(model2, cov.prior = "correlation(gamma, inverse.wishart)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "correlation"));
    checkEquals(length(prior[[1]]@families), d.1 + 1);
    checkEquals(length(prior[[1]]@scales), d.1 + 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 * d.1 + 2 + d.1^2);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "correlation"));
    checkEquals(length(prior[[2]]@families), d.2 + 1);
    checkEquals(length(prior[[2]]@scales), d.2 + 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 * d.2 + 2 + d.2^2);

    checkEquals(prior[[1]]@families[1:d.1], rep(getEnumOrder(familyEnum, "gamma"), d.1));
    checkEquals(prior[[1]]@families[d.1 + 1], getEnumOrder(familyEnum, "inverse.wishart"));
    
    checkEquals(prior[[2]]@families[1:d.2], rep(getEnumOrder(familyEnum, "gamma"), d.2));
    checkEquals(prior[[2]]@families[d.2 + 1], getEnumOrder(familyEnum, "inverse.wishart"));
  }

  prior <- parsePrior(model2, cov.prior = "correlation(gamma(common.scale = 'false'), inverse.wishart)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "correlation"));
    checkEquals(length(prior[[1]]@families), d.1 + 1);
    checkEquals(length(prior[[1]]@scales), d.1 + 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 * d.1 + 2 + d.1^2);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "correlation"));
    checkEquals(length(prior[[2]]@families), d.2 + 1);
    checkEquals(length(prior[[2]]@scales), d.2 + 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 * d.2 + 2 + d.2^2);

    checkEquals(prior[[1]]@families[1:d.1], rep(getEnumOrder(familyEnum, "gamma"), d.1));
    checkEquals(prior[[1]]@families[d.1 + 1], getEnumOrder(familyEnum, "inverse.wishart"));

    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, blme:::defaultCorrelationPosteriorScale),
                            getEnumOrder(commonScaleEnum, "false"));
    checkEquals(prior[[1]]@scales[1:d.1], rep(scaleInt, d.1));
    
    checkEquals(prior[[2]]@families[1:d.2], rep(getEnumOrder(familyEnum, "gamma"), d.2));
    checkEquals(prior[[2]]@families[d.2 + 1], getEnumOrder(familyEnum, "inverse.wishart"));
    checkEquals(prior[[2]]@scales[1:d.2], rep(scaleInt, d.2));
  }
  
  prior <- parsePrior(model2, cov.prior =
                      "g.1 ~ correlation(inverse.gamma, wishart), g.2 ~ spectral(gamma)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "correlation"));
    checkEquals(length(prior[[1]]@families), d.1 + 1);
    checkEquals(length(prior[[1]]@scales), d.1 + 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 * d.1 + 2 + d.1^2);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "spectral"));
    checkEquals(length(prior[[2]]@families), d.2);
    checkEquals(length(prior[[2]]@scales), d.2);
    checkEquals(length(prior[[2]]@hyperparameters), 2 * d.2);
                
    checkEquals(prior[[1]]@families[1:d.1], rep(getEnumOrder(familyEnum, "inverse.gamma"), d.1));
    checkEquals(prior[[1]]@families[d.1 + 1], getEnumOrder(familyEnum, "wishart"));
    
    checkEquals(prior[[2]]@families[1:d.2], rep(getEnumOrder(familyEnum, "gamma"), d.2));
  }
  
  prior <- parsePrior(model2, cov.prior = "g.1 ~ wishart, g.2 ~ inverse.wishart");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + d.1^2);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + d.2^2);

    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
  }

  prior <- parsePrior(model2, cov.prior = "wishart");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + d.1^2);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + d.2^2);

    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "wishart"));
  }
  
  prior <- parsePrior(model2, cov.prior = "wishart(df = expression(factorDimension + 4))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + d.1^2);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + d.2^2);

    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[1]]@hyperparameters[1], d.1 + 4);
    
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@hyperparameters[1], d.2 + 4);
  }

  prior <- parsePrior(model2, cov.prior = "wishart(scale = expression(diag(dataSd^2, factorDimension)))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + d.1^2);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + d.2^2);

    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[1]]@hyperparameters[1:(d.1^2) + 2], as.numeric(diag(1 / var(y), d.1)));
    
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@hyperparameters[1:(d.2^2) + 2], as.numeric(diag(1 / var(y), d.2)));
  }

  prior <- parsePrior(model2, cov.prior = "wishart(common.scale = 'false')");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + d.1^2);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + d.2^2);

    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "wishart"));

    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, blme:::VAR_SCALE_NAME),
                            getEnumOrder(commonScaleEnum, "false"));

    checkEquals(prior[[1]]@scales[1], scaleInt);
    checkEquals(prior[[2]]@scales[1], scaleInt);
  }

  options(warn = -1);
  model3 <- blmer(y ~ x.1 + (1 + x.1 | g.1) + (1 | g.2),
                  control = list(maxIter = 0L),
                  cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);
  d.1 <- 2;
  d.2 <- 1;
  
  prior <- parsePrior(model3, cov.prior = "inverse.gamma, wishart");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2);

    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "inverse.gamma"));
  }
  
  prior <- parsePrior(model3, cov.prior = "wishart(df = expression(factorDimension + 1))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[1]]@hyperparameters[1], 3);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + 1);
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@hyperparameters[1], 2);
  }
  
  prior <- parsePrior(model3, cov.prior = "wishart(scale = 0.01)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[1]]@hyperparameters[3:6], as.numeric(diag(100.0, 2)));
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + 1);
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@hyperparameters[3], 100.0);
  }
  
  prior <- parsePrior(model3, cov.prior = "wishart(scale = expression(diag(0.01, factorDimension)))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[1]]@hyperparameters[3:6], as.numeric(diag(100.0, 2)));
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + 1);
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@hyperparameters[3], 100.0);
  }
  
  prior <- parsePrior(model3, cov.prior = "g.1 ~ wishart");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
  }
  
  prior <- parsePrior(model3, cov.prior = "g.1 ~ wishart(df = 3)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[1]]@hyperparameters[1], 3);
  }
  
  prior <- parsePrior(model3, cov.prior = "g.1 ~ wishart(scale = 0.01)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[1]]@hyperparameters[3:6], as.numeric(diag(100.0, 2)));
  }
  
  prior <- parsePrior(model3, cov.prior = "g.1 ~ wishart(scale = diag(0.01, 2))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[1]]@hyperparameters[3:6], as.numeric(diag(100.0, 2)));
  }
  
  
  prior <- parsePrior(model3, cov.prior = "inverse.wishart(df = expression(factorDimension + 1))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
    checkEquals(prior[[1]]@hyperparameters[1], 3);
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + 1);
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
    checkEquals(prior[[2]]@hyperparameters[1], 2);
  }
  
  prior <- parsePrior(model3, cov.prior = "inverse.wishart(inverse.scale = 0.01)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
    checkEquals(prior[[1]]@hyperparameters[3:6], as.numeric(diag(0.01, 2)));
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + 1);
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
    checkEquals(prior[[2]]@hyperparameters[3], 0.01);
  }
  
  prior <- parsePrior(model3, cov.prior = "inverse.wishart(inverse.scale = expression(diag(0.01, factorDimension)))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
    checkEquals(prior[[1]]@hyperparameters[3:6], as.numeric(diag(0.01, 2)));
    
    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + 1);
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
    checkEquals(prior[[2]]@hyperparameters[3], 0.01);
  }
  
  prior <- parsePrior(model3, cov.prior = "g.1 ~ inverse.wishart");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
  }
  
  prior <- parsePrior(model3, cov.prior = "g.1 ~ inverse.wishart(df = 3)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@hyperparameters[1], 3);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
  }
  
  prior <- parsePrior(model3, cov.prior = "g.1 ~ inverse.wishart(inverse.scale = 0.01)");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@hyperparameters[3:6], as.numeric(diag(0.01, 2)));
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
  }
  
  prior <- parsePrior(model3, cov.prior = "g.1 ~ inverse.wishart(inverse.scale = diag(0.01, 2))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "inverse.wishart"));
    checkEquals(prior[[1]]@hyperparameters[3:6], as.numeric(diag(0.01, 2)));
  }

  prior <- parsePrior(model3, cov.prior = "wishart(data.scale = 'free', scale = 1, df = expression(factorDimension + 3))");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + d.1^2);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    testCovar <- diag(c(var(y), var(y) / var(x.1)), 2);
    checkEquals(prior[[1]]@hyperparameters,
                c(d.1 + 3, determinant(testCovar, logarithm="true")$modulus,
                  as.numeric(solve(testCovar))));

    checkEquals(prior[[2]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[2]]@families), 1);
    checkEquals(length(prior[[2]]@scales), 1);
    checkEquals(length(prior[[2]]@hyperparameters), 2 + d.2^2);
    checkEquals(prior[[2]]@families[1], getEnumOrder(familyEnum, "wishart"));
    checkEquals(prior[[2]]@hyperparameters,
                c(d.2 + 3, log(var(y)), 1 / var(y)));
  }

  prior <- parsePrior(model3, cov.prior = "g.1 ~ wishart(common.scale = 'false')");
  {
    checkEquals(prior[[1]]@type, getEnumOrder(typeEnum, "direct"));
    checkEquals(length(prior[[1]]@families), 1);
    checkEquals(length(prior[[1]]@scales), 1);
    checkEquals(length(prior[[1]]@hyperparameters), 2 + 4);
    checkEquals(prior[[1]]@families[1], getEnumOrder(familyEnum, "wishart"));
    scaleInt <- getScaleInt(getEnumOrder(posteriorScaleEnum, "var"), # sd not making sense for matrix at the moment
                            getEnumOrder(commonScaleEnum, "false"));
    checkEquals(prior[[1]]@scales, scaleInt);
  }
}
