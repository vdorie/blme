cat("\n\nRUnit test cases for bmer::blmer function with priors on the common scale\n\n");

test.bmer.blmer.varPrior <- function()
{
  # options
  generateData <- FALSE;
  cacheOptimizations <- TRUE;
  
  testRoot <- file.path(.path.package(package="blme"), "unitTests");
  
  # need this to convert st matrices in lmer fit to a vector
  # for comparison against an optim call
  source(file.path(testRoot, "alternateImplementation_matrix.R"), TRUE);
  
  if (!cacheOptimizations) {
    source(file.path(testRoot, "alternateImplementation.R"), TRUE);
    source(file.path(testRoot, "alternateImplementation_commonScale.R"), TRUE);
    source(file.path(testRoot, "alternateImplementation_priors.R"), TRUE);
  }
  
  # eight schools
  y <- c(28, 8, -3, 7, -1, 1, 18, 12);
  sigma <- c(15, 10, 16, 11, 9, 11, 10, 18);

  y.z <- (y - mean(y)) / sigma;

  g <- 1:8
  
  eightSchools <- blmer(y.z ~ 1 + (1 | g), var.prior = "point(1)",
                        cov.prior = NULL, fixef.prior = NULL);

  checkEquals(eightSchools@ST[[1]][1], 0);

  source(file.path(testRoot, "lmmData.R"), TRUE);
  
  options(warn = -1);
  testModel <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                     control=list(maxIter=0L),
                     cov.prior = NULL, fixef.prior = NULL, var.prior = "inverse.gamma(2, 1.0)");
  options(warn = 0);

  if (FALSE) {
    lowerBounds <- unlist(sapply(testModel@ST, function(m) { n <- nrow(m); cons1 <- rep(0, n); cons2 <- rep(-Inf, n * (n - 1) / 2); return(c(cons1, cons2)) }))
    upperBounds <- unlist(sapply(testModel@ST, function(m) { n <- nrow(m); return(rep(Inf, n * (n + 1) / 2)) } ))
    
    optimResults <- optim(blme_stMatricesToVector(testModel@ST), blme_getObjectiveFunctionForParameters,
                          lower=lowerBounds, upper=upperBounds,
                          method="L-BFGS-B", model = testModel, control=list(factr=1e-10));
  } else {
    optimResults <- list(par = c(0.723836301225175, 2.17974117787006, -0.377201558806748, 0.79632717532886, 0.741731901625206, 0, -0.402480069027705, -0.82484790952046, 0.693673747070969));
  }

    blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                      cov.prior = NULL, fixef.prior = NULL, var.prior = "inverse.gamma(2, 1.0)");

  checkEquals(blme_stMatricesToVector(blmerFit@ST), optimResults$par, tolerance=1e-5);
}
