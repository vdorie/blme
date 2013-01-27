cat("\n\nRUnit test cases for bmer:::blmer function with priors on the covariance of the modeled coefficients.\n\n");

test.blme.blmer.covarPrior <- function()
{
  # options
  generateData <- FALSE;
  cacheOptimizations <- TRUE;
  
  testRoot <- file.path(.path.package(package="blme"), "unitTests");
  
  source(file.path(testRoot, "lmmData.R"), TRUE);
  
  # need this to convert st matrices in lmer fit to a vector
  # for comparison against an optim call
  source(file.path(testRoot, "alternateImplementation_matrix.R"), TRUE);
  
  if (!cacheOptimizations) {
    source(file.path(testRoot, "alternateImplementation.R"), TRUE);
    source(file.path(testRoot, "alternateImplementation_commonScale.R"), TRUE);
    source(file.path(testRoot, "alternateImplementation_priors.R"), TRUE);
  }

  getLowerBounds <- function(model) {
    unlist(sapply(model@ST, function(m) {
      n <- nrow(m);
      car <- rep(0, n);
      cdr <- rep(-Inf, n * (n - 1) / 2);
      return(c(car, cdr)) }));
  }
  getUpperBounds <- function(model) {
    unlist(sapply(model@ST, function(m) {
      n <- nrow(m);
      return(rep(Inf, n * (n + 1) / 2)); }));
  }

  
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                     control=list(maxIter=0L), var.prior = NULL,
                     cov.prior = NULL, fixef.prior = NULL);
  options(warn = 0);

  if (!cacheOptimizations) {
    lowerBounds <- getLowerBounds(testModel);
    upperBounds <- getUpperBounds(testModel);

    startingParameters <- blme_stMatricesToVector(testModel@ST);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters, model = testModel,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.685799044401953, 2.08564789987389, -0.373384163195289, 0.756409679230464, 0.715198281606338, 0, -0.42926482121343, -0.82061092083679, 0.69045118177096));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                    cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);
    
  checkEquals(blme_stMatricesToVector(blmerFit@ST), optimResults$par, tolerance=1.5e-5);


  cov.prior <- "g.1 ~ gamma(rate = 0.5)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 | g.1), control=list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);

  if (!cacheOptimizations) {
    lowerBounds <- getLowerBounds(testModel);
    upperBounds <- getUpperBounds(testModel);
    
    startingParameters <- blme_stMatricesToVector(testModel@ST);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.617248561865064));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 | g.1),
                    cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);

  checkEquals(blme_stMatricesToVector(blmerFit@ST), optimResults$par, tolerance=1e-5);
  
  

  cov.prior <- "g.1 ~ inverse.gamma(scale = 2.0)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 | g.1), control = list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);
  
  if (!cacheOptimizations) {
    lowerBounds <- getLowerBounds(testModel);
    upperBounds <- getUpperBounds(testModel);
    
    startingParameters <- blme_stMatricesToVector(testModel@ST);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.72329764027013));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 | g.1),
                    cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);

  checkEquals(blme_stMatricesToVector(blmerFit@ST), optimResults$par, tolerance=1e-5);
  
  
  
  
  cov.prior <- "g.1 ~ wishart(scale = 2)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), control = list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);
  
  if (!cacheOptimizations) {
    lowerBounds <- getLowerBounds(testModel);
    upperBounds <- getUpperBounds(testModel);
    
    startingParameters <- blme_stMatricesToVector(testModel@ST);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.826102911329016, 1.46347018293258, -0.0698662045097548));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1),
                    cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);

  checkEquals(blme_stMatricesToVector(blmerFit@ST), optimResults$par, tolerance=1e-6);
  

  cov.prior <- "g.1 ~ inverse.wishart(inverse.scale = 2)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), control=list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);
  
  if (!cacheOptimizations) {
    lowerBounds <- getLowerBounds(testModel);
    upperBounds <- getUpperBounds(testModel);

    startingParameters <- blme_stMatricesToVector(testModel@ST);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.674426918784064, 1.04989476142012, -0.0579534956984132));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1),
                    cov.prior = "g.1 ~ inverse.wishart(inverse.scale = 2)",
                    fixef.prior = NULL, var.prior = NULL);

  checkEquals(blme_stMatricesToVector(blmerFit@ST), optimResults$par, tolerance=1e-5);



  cov.prior <- "g.1 ~ inverse.gamma(scale = 0.5, posterior.scale = 'var', common.scale = FALSE)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 | g.1), control=list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);

  if (!cacheOptimizations) {
    lowerBounds <- getLowerBounds(testModel);
    upperBounds <- getUpperBounds(testModel);
    
    startingParameters <- blme_stMatricesToVector(testModel@ST);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.344982507403144));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 | g.1),
                    cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);

  checkEquals(blme_stMatricesToVector(blmerFit@ST), optimResults$par, tolerance=1e-5);


  cov.prior <- "g.1 ~ gamma(rate = 2, posterior.scale = 'var', common.scale = FALSE)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 | g.1), control=list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);

  if (!cacheOptimizations) {
    lowerBounds <- c(getLowerBounds(testModel), 0.0);
    upperBounds <- c(getUpperBounds(testModel), Inf);
    
    startingParameters <- c(blme_stMatricesToVector(testModel@ST), 1.0);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.344630660995718, 2.10611396651749));
  }

  blmerFit <- blmer(y ~ x.1 + x.2 + (1 | g.1),
                    cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);

  blmerResults <- c(blme_stMatricesToVector(blmerFit@ST), blmerFit@deviance[["sigmaREML"]]);
  checkEquals(blmerResults, optimResults$par, tolerance = 1e-6);


  cov.prior <- "g.1 ~ inverse.gamma(scale = 0.5, posterior.scale = 'sd', common.scale = FALSE)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 | g.1), control=list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);

  if (!cacheOptimizations) {
    lowerBounds <- c(getLowerBounds(testModel), 0.0);
    upperBounds <- c(getUpperBounds(testModel), Inf);
    
    startingParameters <- c(blme_stMatricesToVector(testModel@ST), 1.0);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.384183905569239, 2.09594620590208));
  }

  blmerFit <- blmer(y ~ x.1 + x.2 + (1 | g.1),
                    cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);

  blmerResults <- c(blme_stMatricesToVector(blmerFit@ST), blmerFit@deviance[["sigmaREML"]]);
  checkEquals(blmerResults, optimResults$par, tolerance = 1e-6);


  cov.prior <- "g.1 ~ gamma(rate = 2, posterior.scale = 'sd', common.scale = FALSE)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 | g.1), control=list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);

  if (!cacheOptimizations) {
    lowerBounds <- c(getLowerBounds(testModel), 0.0);
    upperBounds <- c(getUpperBounds(testModel), Inf);
    
    startingParameters <- c(blme_stMatricesToVector(testModel@ST), 1.0);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.409482705974117, 2.09078910076072));
  }

  blmerFit <- blmer(y ~ x.1 + x.2 + (1 | g.1),
                    cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);

  blmerResults <- c(blme_stMatricesToVector(blmerFit@ST), blmerFit@deviance[["sigmaREML"]]);
  checkEquals(blmerResults, optimResults$par, tolerance = 1e-6);


  cov.prior <- "g.1 ~ wishart(scale = 2, common.scale = FALSE)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), control = list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);
  
  if (!cacheOptimizations) {
    lowerBounds <- c(getLowerBounds(testModel), 0.0);
    upperBounds <- c(getUpperBounds(testModel), Inf);
    
    startingParameters <- c(blme_stMatricesToVector(testModel@ST), 1.0);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.792440504958109, 1.37292598078286, -0.0685088488977944, 1.26420602775867));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1),
                    cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);

  blmerResults <- c(blme_stMatricesToVector(blmerFit@ST), blmerFit@deviance[["sigmaREML"]]);
  checkEquals(blmerResults, optimResults$par, tolerance=1e-6);


  cov.prior <- "g.1 ~ inverse.wishart(inverse.scale = 2, common.scale = FALSE)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1), control = list(maxIter=0L),
                     cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);
  
  if (!cacheOptimizations) {
    lowerBounds <- getLowerBounds(testModel);
    upperBounds <- getUpperBounds(testModel);
    
    startingParameters <- blme_stMatricesToVector(testModel@ST);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.59154080359443, 1.03156958901312, -0.0949257670583755));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1),
                    cov.prior = cov.prior, fixef.prior = NULL, var.prior = NULL);

  blmerResults <- blme_stMatricesToVector(blmerFit@ST);
  checkEquals(blmerResults, optimResults$par, tolerance=1e-5);
}
