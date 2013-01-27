getDirectDefault <- function(regression, callingEnvironment, factorNumber, family)
{
  result <- list(family = family);

  fillDefaults <- function(x, ...) x;
  if (result$family == GAMMA_FAMILY_NAME) {
    fillDefaults <- getDirectGammaDefaults;
  } else if (result$family == INVGAMMA_FAMILY_NAME) {
    fillDefaults <- getDirectInverseGammaDefaults;
  } else if (result$family == WISHART_FAMILY_NAME) {
    fillDefaults <- getDirectWishartDefaults;
  } else if (result$family == INVWISHART_FAMILY_NAME) {
    fillDefaults <- getDirectInverseWishartDefaults;
  }

  prior <- fillDefaults(regression, callingEnvironment, factorNumber, result);
  if (factorNumber > 0) return(prior);
  
  return(list(prior = prior, fillDefaults = fillDefaults));
}

parseDirectSpecification <- function(regression, callingEnvironment, factorNumber, family,
                                     specification)
{
  if (factorNumber <= 0) {
    prior <- getDirectDefault(regression, callingEnvironment, factorNumber, family);
  } else {
    factorDimension <- nrow(regression@ST[[factorNumber]]);

    if (factorDimension > 1 && any(family == UNIVARIATE_FAMILY_NAMES)) {
      stop(family, " applied to multivariate factor '",
              getFactorNameForNumber(regression, factorNumber));
    }
    # because the general machinery is to support multiple coordinate priors
    # per factor, we store the prior for the factor in a list of size 1
    prior <- list(getDirectDefault(regression, callingEnvironment, factorNumber, family));
  }
  
  option <- trim(specification);
  if (nchar(option) > 0) {
    prior <- parseDirectOption(regression, callingEnvironment, factorNumber, family, option);
  }

  # Turn the result into something that resembles the priors created by
  # the correlation and spectral priors, namely lists of lists, with
  # possible defaults for individual coordinates.
  result <- NULL;
  if (factorNumber == 0) {
    result <- prior;
    result$type <- DIRECT_TYPE_NAME;
  } else {
    result <- list(type = DIRECT_TYPE_NAME, prior = prior, coordinateDefault = NULL);
    
    return(result);
  }
  
  return (result);
}


parseDirectOption <- function(regression, callingEnvironment, factorNumber, family,
                              options)
{
  optionsList <- eval(parse(text=paste("list(", options, ")", sep="")),
                      envir=callingEnvironment);
  
  prior <- completeDirectPriorFromOptions(regression, callingEnvironment,
                                          family, optionsList, factorNumber);

  if (factorNumber <= 0) {
    return(prior);
  } else {
    return(list(prior));
  }
}

completeDirectPriorFromOptions <- function(regression, callingEnvironment,
                                           family, optionsList, factorNumber)
{
  if (is.null(names(optionsList))) {
    namedOptions   <- list();
    unnamedOptions <- optionsList;
  } else {
    namedOptions   <- optionsList[names(optionsList) != ""];
    unnamedOptions <- optionsList[names(optionsList) == ""];
  }

  # references
  namedOptionsRef   <- list(env = sys.frame(sys.nframe()), var = "namedOptions");
  unnamedOptionsRef <- list(env = sys.frame(sys.nframe()), var = "unnamedOptions");
  
  prior <- list(family=family);
  fillDefaults <- function(x, ...) x;

  if (family == GAMMA_FAMILY_NAME) {
    prior$shape <- getPriorOption(SHAPE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    rate <- getPriorOption(RATE_HYPERPARAMETER_NAME,  namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$commonScale    <- getPriorOption(COMMON_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    # if the user sets a rate but not a data scale, force data scale to absolute
    if (!is.null(rate) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$rate <- rate;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getDirectGammaDefaults;
  } else if (family == INVGAMMA_FAMILY_NAME) {
    prior$shape <- getPriorOption(SHAPE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    scale <- getPriorOption(SCALE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$commonScale    <- getPriorOption(COMMON_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    if (!is.null(scale) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$scale <- scale;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getDirectInverseGammaDefaults;
  } else if (family == WISHART_FAMILY_NAME) {
    prior$degreesOfFreedom <-
      getPriorOption(DEGREES_OF_FREEDOM_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    scale <- getPriorOption(SCALE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$commonScale    <- getPriorOption(COMMON_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    if (!is.null(scale) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$scale <- scale;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getDirectWishartDefaults;
  } else if (family == INVWISHART_FAMILY_NAME) {
    prior$degreesOfFreedom <-
      getPriorOption(DEGREES_OF_FREEDOM_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    inverseScale <- getPriorOption(INVERSE_SCALE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$commonScale    <- getPriorOption(COMMON_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    if (!is.null(inverseScale) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$inverseScale <- inverseScale;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getDirectInverseWishartDefaults;
  }

  if (length(namedOptions) > 0) {
    warning("Unrecognized prior option(s) for ", family, " family: ",
            toString(names(namedOptions)), ".");
  }
  if (length(unnamedOptions) > 0) {
    warning("Extra option(s) for ", family, " family: ", toString(unnamedOptions), ".");
  }

  prior <- fillDefaults(regression, callingEnvironment, factorNumber, prior);
  if (factorNumber > 0) return(prior);
  return(list(prior = prior, fillDefaults = fillDefaults));
}

getDirectGammaDefaults <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (is.null(prior$shape))          prior$shape <- defaultDirectGammaShape;
  if (is.null(prior$posteriorScale)) prior$posteriorScale <- defaultDirectPosteriorScale;
  if (is.null(prior$commonScale))    prior$commonScale <- defaultDirectCommonScale;
  if (is.null(prior$rate))           prior$rate <- defaultDirectGammaRate;
  if (is.null(prior$dataScale))      prior$dataScale <- defaultDirectDataScale;
  
  return(evaluateGammaExpressions(regression, callingEnvironment, factorNumber, prior));
}

getDirectInverseGammaDefaults <- function(regression, callingEnvironment, factorNumber, prior) {
  if (is.null(prior$shape))          prior$shape <- defaultDirectInverseGammaShape;
  if (is.null(prior$posteriorScale)) prior$posteriorScale <- defaultDirectPosteriorScale;
  if (is.null(prior$commonScale))    prior$commonScale <- defaultDirectCommonScale;
  if (is.null(prior$scale))          prior$scale <- defaultDirectInverseGammaScale;
  if (is.null(prior$dataScale))      prior$dataScale <- defaultDirectDataScale;

  return(evaluateInverseGammaExpressions(regression, callingEnvironment, factorNumber, prior));
}

getDirectWishartDefaults <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (is.null(prior$degreesOfFreedom))
    prior$degreesOfFreedom <- defaultDirectWishartDegreesOfFreedom;
  if (is.null(prior$scale))
    prior$scale <- defaultDirectWishartScale;
  prior$posteriorScale <- VAR_SCALE_NAME; # sqrt of matrices not currently supported
  if (is.null(prior$commonScale))
    prior$commonScale <- defaultDirectCommonScale;
  if (is.null(prior$dataScale))
    prior$dataScale <- defaultDirectDataScale;
  
  return(evaluateWishartExpressions(regression, callingEnvironment, factorNumber, prior));
}

getDirectInverseWishartDefaults <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (is.null(prior$degreesOfFreedom))
    prior$degreesOfFreedom <- defaultDirectInverseWishartDegreesOfFreedom;
  if (is.null(prior$inverseScale))
    prior$inverseScale <- defaultDirectInverseWishartInverseScale;
  prior$posteriorScale <- VAR_SCALE_NAME; # sqrt of matrices not currently supported
  if (is.null(prior$commonScale))
    prior$commonScale <- defaultDirectCommonScale;
  if (is.null(prior$dataScale))
    prior$dataScale <- defaultDirectDataScale;

  return(evaluateInverseWishartExpressions(regression, callingEnvironment, factorNumber, prior));
}
