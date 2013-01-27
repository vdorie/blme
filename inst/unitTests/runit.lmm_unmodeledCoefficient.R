cat("\n\nRUnit test cases for bmer::blmer function with priors on the unmodeled coefficients\n\n");

test.bmer.blmer.fixefPrior <- function()
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
                     control=list(maxIter=0L),
                     cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);
  options(warn = 0);

  if (!cacheOptimizations) {
    lowerBounds <- getLowerBounds(testModel);
    upperBounds <- getUpperBounds(testModel);

    startingParameters <- blme_stMatricesToVector(testModel@ST);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <-
      list(par = c(0.685799044401953, 2.08564789987389, -0.373384163195289, 0.756409679230464, 0.715198281606338, 0, -0.42926482121343, -0.82061092083679, 0.69045118177096));
  }
    
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                    cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);

  blmerResults <- blme_stMatricesToVector(blmerFit@ST);
  checkEquals(blmerResults, optimResults$par, tolerance=1.5e-5);


  fixef.prior <- "normal(sd = 1, common.scale = TRUE)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                     control=list(maxIter=0L),
                     cov.prior = NULL, var.prior = NULL, fixef.prior = fixef.prior);
  options(warn = 0);

  if (!cacheOptimizations) {
    lowerBounds <- getLowerBounds(testModel);
    upperBounds <- getUpperBounds(testModel);
    
    startingParameters <- blme_stMatricesToVector(testModel@ST);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.645491684136578, 1.94805334962143, -0.42573815389605,
                           5.24061970057903, 0.824637022004438, 0.54266484181711,
                           0.278547275225769, 0.768154464529237, 1.42827591442321));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                    cov.prior = NULL, var.prior = NULL, fixef.prior = fixef.prior);

  blmerResults <- blme_stMatricesToVector(blmerFit@ST);
  checkEquals(blmerResults, optimResults$par, tolerance=1e-5);


  fixef.prior <- "normal(sd = 1, common.scale = FALSE)";
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                     control=list(maxIter=0L),
                     cov.prior = NULL, var.prior = NULL, fixef.prior = fixef.prior);
  options(warn = 0);

  if (!cacheOptimizations) {
    lowerBounds <- c(getLowerBounds(testModel), 0.0);
    upperBounds <- c(getUpperBounds(testModel), Inf);
    
    startingParameters <- c(blme_stMatricesToVector(testModel@ST), 1.0);
    optimResults <- optim(startingParameters, blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <-
      list(par = c(0.638289608865955, 1.93169351638693, -0.420471278094413, 5.04396706905165,
             0.831562587802851, 0.496972740884312, 0.291250499210988, 0.768201439493807,
             1.41639282835109, 0.923663192383467));
  }
  
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                    cov.prior = NULL, var.prior = NULL, fixef.prior = fixef.prior);

  blmerResults <- c(blme_stMatricesToVector(blmerFit@ST), blmerFit@deviance[["sigmaREML"]]);
  checkEquals(blmerResults, optimResults$par, tolerance=1e-5);
  
  
  # check even that these run without error
  options(warn = 2);
  ignored <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
                   cov.prior = NULL, var.prior = NULL,
                   fixef.prior = "normal");
  ignored <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
                   cov.prior = NULL, var.prior = NULL,
                   fixef.prior = "normal(cov = diag(0.5, 2), data.scale='absolute')");
  options(warn = 0);
}
