## hack this as on "common scale" is really inconvenient here
getResidPriorDFAdjustment <- function(residPrior)
{
  if (is(residPrior, "bmerGammaDist")) {
    return(-(residPrior@shape - 1.0) * if (residPrior@posteriorScale == 'sd') 1 else 2);
  } else if (is(residPrior, "bmerInvGammaDist")) {
    return( (residPrior@shape + 1.0) * if (residPrior@posteriorScale == 'sd') 1 else 2);
  }
  return(0);
}
    
createBlmerControl <- function(reTrms, priors)
{
  df <- 0; ## adjustment to polynomial (sigma.sq)^{-df/2}
  constant <- 0; ## normalizing constants and the like. On deviance (-2 log) scale

  numFactors <- length(priors$covPriors);
 
  df <- df + getDFAdjustment(priors$fixefPrior) + getResidPriorDFAdjustment(priors$residPrior);
  constant <- constant + getConstantTerm(priors$fixefPrior) + getConstantTerm(priors$residPrior);

  for (i in 1:numFactors) {
    df <- df + getDFAdjustment(priors$covPrior[[i]]);
    constant <- constant + getConstantTerm(priors$covPrior[[i]]);
  }

  sigmaOptimizationType <- getSigmaOptimizationType(priors);

  numCovParameters <- sum(sapply(reTrms$cnms, function(cnm) { d <- length(cnm); d * (d + 1) / 2; }));

  ranefStructure <- list(numRanefPerFactor = diff(reTrms$Gp),
                         numCoefPerFactor = as.integer(sapply(reTrms$cnms, length)),
                         numFactors = length(reTrms$cnms));
  ranefStructure$numGroupsPerFactor <- as.integer(ranefStructure$numRanefPerFactor / ranefStructure$numCoefPerFactor + 0.5);
  
  return(list(df = df, constant = constant,
              sigmaOptimizationType = sigmaOptimizationType,
              numCovParameters = numCovParameters,
              ranefStructure = ranefStructure));
}

SIGMA_OPTIM_NUMERIC      <- "numeric";
SIGMA_OPTIM_POINT        <- "point";
SIGMA_OPTIM_SQ_LINEAR    <- "sigma.sq.linear";
SIGMA_OPTIM_SQ_QUADRATIC <- "sigma.sq.quadratic";
SIGMA_OPTIM_QUADRATIC    <- "sigma.quadratic";
## determines how to optimize over sigma (and maybe other stuff in future)
## possible values are:
##  sigmaOptimization ==
##    SIGMA_OPTIM_NUMERIC      - brute force by adding to numeric optimizer
##    SIGMA_OPTIM_POINT        - sigma is fixed to a particular value
##    SIGMA_OPTIM_SQ_LINEAR    - sigma.sq.hat is root to linear equation
##    SIGMA_OPTIM_SQ_QUADRATIC - sigma.sq.hat is root to quadratic equation
##    SIGMA_OPTIM_QUADRATIC    - sigma.hat is root to quadratic equation
getSigmaOptimizationType <- function(priors)
{
  fixefPrior <- priors$fixefPrior;
  covPriors  <- priors$covPriors;
  residPrior <- priors$residPrior;

  if (is(residPrior, "bmerPointDist")) return(SIGMA_OPTIM_POINT);
  
  if (is(fixefPrior, "bmerNormalDist") && fixefPrior@commonScale == FALSE)
    return(SIGMA_OPTIM_NUMERIC);

  
  exponentialTerms <- c();
  for (i in 1:length(covPriors)) {
    covPrior.i <- covPriors[[i]];
    exponentialTerm <- getExponentialSigmaPower(covPrior.i);
    if (exponentialTerm != 0) exponentialTerms <- union(exponentialTerms, exponentialTerm);
  }
  exponentialTerm <- getExponentialSigmaPower(residPrior);
  if (exponentialTerm != 0) exponentialTerms <- union(exponentialTerms, exponentialTerm);

  ## exp(-0.5 * sigma^-2 * stuff) always happens, so other terms are "extra"
  extraExponentialTerms <- setdiff(exponentialTerms, -2);
  
  if (length(extraExponentialTerms) == 0) return(SIGMA_OPTIM_SQ_LINEAR);
  
  if (length(extraExponentialTerms) > 1 || !(extraExponentialTerms %in% c(-1, 2)))
    return(SIGMA_OPTIM_NUMERIC);

  if (extraExponentialTerms == -1) return(SIGMA_OPTIM_QUADRATIC);
  if (extraExponentialTerms ==  2) return(SIGMA_OPTIM_SQ_QUADRATIC);

  ## should be unreachable
  return(SIGMA_OPTIM_NUMERIC);
}
