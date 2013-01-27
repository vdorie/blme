###
# The following set of functions flesh out any priors passed to them by
# creating and populating an expression and using that to evaluate any
# expressions acting as placeholders in the priors. At present, expresions
# are only evaluated when the prior is associated with a specific factor.
##
evaluateGammaExpressions <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (factorNumber <= 0) return(prior);
  
  env <- new.env(parent = callingEnvironment);
  env$dataSd <- ifelse(isLinearMixedModel(regression), sd(regression@y), 1);
  env$inputSd <- getInputSdForFactor(regression, factorNumber);
  env$sdRatio <- env$dataSd / env$inputSd;
  
  loadConstants(env);

  # the checks would normally be in the validate function, but we need to know this before
  # we can evaluate some subsequent expressions
  errorPrefix <- paste("Error applying prior to factor '",
                       getFactorNameForNumber(regression, factorNumber),
                       "': ", sep='');
  if (is.expression(prior$posteriorScale))
    prior$posteriorScale <- eval(prior$posteriorScale, env);
  if (prior$posteriorScale != SD_SCALE_NAME &&
      prior$posteriorScale != VAR_SCALE_NAME) {
    stop(errorPrefix,  "posterior scale for univariate prior must be '",
         SD_SCALE_NAME, "' or '", VAR_SCALE_NAME, "', was '",
         prior$posteriorScale, "'.");
  }
  env$posteriorScale <- prior$posteriorScale;
  
  if (is.expression(prior$commonScale))
    prior$commonScale <- eval(prior$commonScale, env);
  if (is.logical(prior$commonScale))
    prior$commonScale <- ifelse(prior$commonScale, COMMON_SCALE_TRUE_NAME, COMMON_SCALE_FALSE_NAME);
  
  if (is.expression(prior$shape))
    prior$shape <- eval(prior$shape, env);
  env$shape <- prior$shape;
  
  if (is.expression(prior$rate))
    prior$rate <- eval(prior$rate, env);

  return(prior);
}

evaluateInverseGammaExpressions <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (factorNumber <= 0) return(prior);
  env <- new.env(parent = callingEnvironment);
  env$dataSd <- ifelse(isLinearMixedModel(regression), sd(regression@y), 1);
  env$inputSd <- getInputSdForFactor(regression, factorNumber);
  env$sdRatio <- env$dataSd / env$inputSd;
  
  loadConstants(env);
  
  errorPrefix <- paste("Error applying prior to factor '",
                       getFactorNameForNumber(regression, factorNumber),
                       "': ", sep='');
  if (is.expression(prior$posteriorScale))
    prior$posteriorScale <- eval(prior$posteriorScale, env);
  if (prior$posteriorScale != SD_SCALE_NAME &&
      prior$posteriorScale != VAR_SCALE_NAME) {
    stop(errorPrefix,  "posterior scale for univariate prior must be '",
         SD_SCALE_NAME, "' or '", VAR_SCALE_NAME, "', was '",
         prior$posteriorScale, "'.");
  }
  env$posteriorScale <- prior$posteriorScale;

  if (is.expression(prior$commonScale))
    prior$commonScale <- eval(prior$commonScale, env);
  if (is.logical(prior$commonScale))
    prior$commonScale <- ifelse(prior$commonScale, COMMON_SCALE_TRUE_NAME, COMMON_SCALE_FALSE_NAME);
    
  if (is.expression(prior$shape))
    prior$shape <- eval(prior$shape, env);
  env$shape <- prior$shape;
  
  if (is.expression(prior$scale))
    prior$scale <- eval(prior$scale, env);

  return(prior);
}

evaluateWishartExpressions <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (factorNumber <= 0) return(prior);
  
  env <- new.env(parent = callingEnvironment);
  env$dataSd <- ifelse(isLinearMixedModel(regression), sd(regression@y), 1);
  env$inputSd <- getInputSdForFactor(regression, factorNumber);
  env$sdRatio <- env$dataSd / env$inputSd;
  env$factorDimension <- nrow(regression@ST[[factorNumber]]);

  if (is.expression(prior$posteriorScale))
    prior$posteriorScale <- eval(prior$posteriorScale, env);
  if (is.expression(prior$commonScale))
    prior$commonScale <- eval(prior$commonScale, env);
  if (is.logical(prior$commonScale))
    prior$commonScale <- ifelse(prior$commonScale, COMMON_SCALE_TRUE_NAME, COMMON_SCALE_FALSE_NAME);
  
  loadConstants(env);
  
  if (is.expression(prior$degreesOfFreedom))
    prior$degreesOfFreedom <- eval(prior$degreesOfFreedom, env);
  env$degreesOfFreedom <- prior$degreesOfFreedom;
  
  if (is.expression(prior$scale))
    prior$scale <- eval(prior$scale, env);
  
  return(prior);
}

evaluateInverseWishartExpressions <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (factorNumber <= 0) return(prior);
  
  env <- new.env(parent = callingEnvironment);
  env$dataSd <- ifelse(isLinearMixedModel(regression), sd(regression@y), 1);
  env$inputSd <- getInputSdForFactor(regression, factorNumber);
  env$sdRatio <- env$dataSd / env$inputSd;
  env$factorDimension <- nrow(regression@ST[[factorNumber]]);

  if (is.expression(prior$posteriorScale))
    prior$posteriorScale <- eval(prior$posteriorScale, env);
  if (is.expression(prior$commonScale))
    prior$commonScale <- eval(prior$commonScale, env);
  if (is.logical(prior$commonScale))
    prior$commonScale <- ifelse(prior$commonScale, COMMON_SCALE_TRUE_NAME, COMMON_SCALE_FALSE_NAME);
  
  loadConstants(env);

  if (is.expression(prior$degreesOfFreedom))
    prior$degreesOfFreedom <- eval(prior$degreesOfFreedom, env);
  env$degreesOfFreedom <- prior$degreesOfFreedom;
  
  if (is.expression(prior$inverseScale))
    prior$inverseScale <- eval(prior$inverseScale, env);
  
  return(prior);
}
