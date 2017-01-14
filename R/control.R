## "control" refers to how optimization should proceed, i.e. which parameters
## need numeric and which can be profiled out

## hack this as on "common scale" is really inconvenient here
getResidPriorDFAdjustment <- function(residPrior)
{
  if (is(residPrior, "bmerGammaDist")) {
    return(-(residPrior@shape - 1.0) * if (residPrior@posteriorScale == 'sd') 1 else 2)
  } else if (is(residPrior, "bmerInvGammaDist")) {
    return( (residPrior@shape + 1.0) * if (residPrior@posteriorScale == 'sd') 1 else 2)
  }
  0
}

getThetaLowerBoundsForDimension <- function(d) {
  if (d == 1) return(0)
  c(0, rep(-Inf, d - 1), getThetaLowerBoundsForDimension(d - 1))
}

## TODO: this should eventually not assume the ranef structure but instead
## suck it from Lind and theta, if possible
getRanefStructure <- function(pred, resp, reTrms) {
  ranefStructure <- list(numCovParameters = sum(sapply(reTrms$cnms, function(cnm) { d <- length(cnm); d * (d + 1) / 2; })),
                         numRanefPerFactor = diff(reTrms$Gp),
                         numCoefPerFactor = as.integer(sapply(reTrms$cnms, length)),
                         numFactors = length(reTrms$cnms))
  ranefStructure$numGroupsPerFactor <- as.integer(ranefStructure$numRanefPerFactor / ranefStructure$numCoefPerFactor + 0.5)
  ranefStructure$lower <- as.numeric(unlist(sapply(ranefStructure$numCoefPerFactor, getThetaLowerBoundsForDimension)))
  
  ranefStructure
}

createBlmerControl <- function(pred, resp, priors)
{
  df <- 0 ## adjustment to polynomial (sigma.sq)^{-df/2}
  constant <- 0 ## normalizing constants and the like. On deviance (-2 log) scale

  numFactors <- length(priors$covPriors)
 
  df <- df + getDFAdjustment(priors$fixefPrior) + getResidPriorDFAdjustment(priors$residPrior)
  constant <- constant + getConstantTerm(priors$fixefPrior) + getConstantTerm(priors$residPrior)

  for (i in seq_len(numFactors)) {
    df <- df + getDFAdjustment(priors$covPrior[[i]])
    constant <- constant + getConstantTerm(priors$covPrior[[i]])
  }

  fixefOptimizationType <- getFixefOptimizationType(pred, resp, priors)
  sigmaOptimizationType <- getSigmaOptimizationType(resp, priors)
  
  namedList(df, constant, fixefOptimizationType, sigmaOptimizationType)
}

FIXEF_OPTIM_NA      <- "na"      ## no fixefs in model
FIXEF_OPTIM_NUMERIC <- "numeric" ## brute force by adding to numeric optimizer
FIXEF_OPTIM_LINEAR  <- "linear"  ## mle found by root of linear equation. or, don't worry about it
getFixefOptimizationType <- function(pred, resp, priors)
{
  if (length(pred$X) == 0) return(FIXEF_OPTIM_NA)
  
  if (!is(resp, "lmerResp")) return(FIXEF_OPTIM_NUMERIC)
  
  fixefPrior <- priors$fixefPrior
  
  if (is(fixefPrior, "bmerTDist") || is(fixefPrior, "bmerHorseshoeDist")) return(FIXEF_OPTIM_NUMERIC)

  FIXEF_OPTIM_LINEAR
}

## determines how to optimize over sigma
## possible values are:
SIGMA_OPTIM_NA           <- "na"                 ## doesn't apply
SIGMA_OPTIM_NUMERIC      <- "numeric"            ## brute force by adding to numeric optimizer
SIGMA_OPTIM_POINT        <- "point"              ## sigma is fixed to a particular value
SIGMA_OPTIM_SQ_LINEAR    <- "sigma.sq.linear"    ## sigma.sq.hat is root to linear equation
SIGMA_OPTIM_SQ_QUADRATIC <- "sigma.sq.quadratic" ## sigma.sq.hat is root to quadratic equation
SIGMA_OPTIM_QUADRATIC    <- "sigma.quadratic"    ## sigma.hat is root to quadratic equation

getSigmaOptimizationType <- function(resp, priors)
{
  if (!is(resp, "lmerResp")) return(SIGMA_OPTIM_NA)
  
  fixefPrior <- priors$fixefPrior
  covPriors  <- priors$covPriors
  residPrior <- priors$residPrior

  if (is(residPrior, "bmerPointDist")) return(SIGMA_OPTIM_POINT)
  
  if (is(fixefPrior, "bmerNormalDist") && fixefPrior@commonScale == FALSE)
    return(SIGMA_OPTIM_NUMERIC)
  if ((is(fixefPrior, "bmerTDist") || is(fixefPrior, "bmerHorseshoeDist")) && fixefPrior@commonScale == TRUE)
    return(SIGMA_OPTIM_NUMERIC)
  
  exponentialTerms <- c()
  for (i in seq_along(covPriors)) {
    covPrior.i <- covPriors[[i]]

    if (is(covPrior.i, "bmerCustomDist") && covPrior.i@commonScale == FALSE) return(SIGMA_OPTIM_NUMERIC)
    
    exponentialTerm <- getExponentialSigmaPower(covPrior.i)
    if (exponentialTerm != 0) exponentialTerms <- union(exponentialTerms, exponentialTerm)
  }
  exponentialTerm <- getExponentialSigmaPower(residPrior)
  if (exponentialTerm != 0) exponentialTerms <- union(exponentialTerms, exponentialTerm)

  ## exp(-0.5 * sigma^-2 * stuff) always happens, so other terms are "extra"
  extraExponentialTerms <- setdiff(exponentialTerms, -2)
  
  if (length(extraExponentialTerms) == 0) return(SIGMA_OPTIM_SQ_LINEAR)
  
  if (length(extraExponentialTerms) > 1 || !(extraExponentialTerms %in% c(-1, 2)))
    return(SIGMA_OPTIM_NUMERIC)

  if (extraExponentialTerms == -1) return(SIGMA_OPTIM_QUADRATIC)
  if (extraExponentialTerms ==  2) return(SIGMA_OPTIM_SQ_QUADRATIC)

  ## should be unreachable
  SIGMA_OPTIM_NUMERIC
}
