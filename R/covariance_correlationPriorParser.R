getCorrelationCoordinateDefault <- function(regression, callingEnvironment, factorNumber)
{
  result <- list(family = defaultCorrelationCoordinateFamily);

  fillDefaults <- function(x, ...) x;
  if (result$family == GAMMA_FAMILY_NAME) {
    fillDefaults <- getCorrelationGammaDefaults;
  } else if (result$family == INVGAMMA_FAMILY_NAME) {
    fillDefaults <- getCorrelationInverseGammaDefaults;
  }

  prior <- fillDefaults(regression, callingEnvironment, factorNumber, result);
  if (factorNumber > 0) return(prior);

  return(list(prior = prior, fillDefaults = fillDefaults));
}

getCorrelationMatrixDefault <- function(regression, callingEnvironment, factorNumber)
{
  result <- list(family = defaultCorrelationMatrixFamily);

  fillDefaults <- function(x, ...) x;
  if (result$family == WISHART_FAMILY_NAME) {
    fillDefaults <- getCorrelationWishartDefaults;
  } else if (result$family == INVWISHART_FAMILY_NAME) {
    fillDefaults <- getCorrelationInverseWishartDefaults;
  }

  prior <- fillDefaults(regression, callingEnvironment, factorNumber, result);
  if (factorNumber > 0) return(prior);

  return(list(prior = prior, fillDefaults = fillDefaults));
}

parseCorrelationSpecification <- function(regression, callingEnvironment,
                                          factorNumber, specification)
{
  if (factorNumber == 0) {
    coordinatePriors <- list();
  } else {
    factorDimension  <- nrow(regression@ST[[factorNumber]]);
    coordinatePriors <- vector("list", factorDimension);
  }
  
  coordinateDefault <- getCorrelationCoordinateDefault(regression, callingEnvironment, factorNumber);
  matrixPrior       <- getCorrelationMatrixDefault(regression, callingEnvironment, factorNumber);
  
  option <- trim(specification);
  while (nchar(option) > 0) {
    # should go through the options one at a time, identifying any defaults
    # and installing the right priors in the right spots. pops off the
    # head of the options comma-separated list
    parseResult <- parseCorrelationOption(regression, callingEnvironment, factorNumber,
                                          coordinatePriors, matrixPrior,
                                          coordinateDefault,
                                          option);
    
    coordinatePriors  <- parseResult$coordinatePriors;
    matrixPrior       <- parseResult$matrixPrior;
    coordinateDefault <- parseResult$coordinateDefault;
    
    option            <- parseResult$option;
  }

  result <- list(type = CORRELATION_TYPE_NAME);
  if (factorNumber == 0) {
    result <- list(type = CORRELATION_TYPE_NAME);
    result$coordinate <- coordinateDefault;
    result$matrix     <- matrixPrior;
  } else {
    priors <- coordinatePriors;
    priors[[factorDimension + 1]] <- matrixPrior;
    
    result$prior = priors;
    result$coordinateDefault = coordinateDefault;
  }

  return (result);
}

parseCorrelationOption <- function(regression, callingEnvironment, factorNumber,
                                   coordinatePriors, matrixPrior, coordinateDefault,
                                   option)
{
  # split the first coordinate prior off from the rest of the comma separated list
  commaSplit <- splitListOnComma(option);
  if (length(commaSplit) == 1) {
    remainder <- "";
  } else {
    option    <- commaSplit[1];
    remainder <- trim(commaSplit[2]);
  }

  result <- list(coordinatePriors =  coordinatePriors,  matrixPrior = matrixPrior,
                 coordinateDefault = coordinateDefault, option      = remainder);

  # Does not support "FLAT", unless someone can convince me that it should.
  supportedFamilies <- flattenStrings(ALL_FAMILY_NAMES, sep="|"); # for regex
  
  isDefault <- FALSE;
  coordinateNumber <- 0;
  # figure out the coordinate name (or none if default), parametric family, and the options
  # force default status if prior is a global default (factorNumber == 0)
  if (isDefaultSpecification(option) || factorNumber == 0) {
    splitFamilyFromOptionsRegex <- paste("^(", supportedFamilies, ")",
                                         "(?:\\s*\\((.*)\\))?$", sep="");
  
    # we use regex with \1 as a special substitution char to split the results
    parse <- subSplit(splitFamilyFromOptionsRegex,
                      "\\1\1\\2",
                      option, perl=TRUE);
    familyName <- parse[1];
    options    <- parse[2];
    rm(parse);

    isDefault <- TRUE;
  } else {
    specificationValidatorRegex <-
      paste("^\\s*\\S+\\s*~\\s*(?:", supportedFamilies, ")", # a ~ gamma
            "\\s*(?:\\(.*\\))?$",                            # a ~ gamma(this = that, what = when)
            sep="");
    
    isValidSpecification <-
      grepl(specificationValidatorRegex, option, perl=TRUE);
    if (!isValidSpecification) stop("Unable to parse '", option, "'.");
    
    splitFactorFamilyAndOptionsRegex <- 
      paste("^\\s*(\\S+)\\s*~\\s*(", supportedFamilies, ")", # same as above, w/ capture groups
            "\\s*(?:\\((.*)\\))?$", sep="");

    parse <- subSplit(splitFactorFamilyAndOptionsRegex,
                      "\\1\1\\2\1\\3",
                      option, perl=TRUE);
    
    coordinateName <- parse[1];
    familyName     <- parse[2];
    options        <- parse[3];
    rm(parse);

    if (any(familyName == UNIVARIATE_FAMILY_NAMES)) {
      coordinateNumber <- tryCatch(as.integer(coordinateName), warning = function(e) NA);
      if (is.na(coordinateNumber)) {
        coordinateNumber <- getCoordinateNumberForName(regression, factorNumber, coordinateName);
        if (coordinateNumber == 0) {
          stop("Unable to identify coordinate name in '", option, "'.");
        }
      }
      if (coordinateNumber < 1 || coordinateNumber > nrow(regression@ST[[factorNumber]]))
        stop("Coordinate number ", coordinateNumber, " out of range.");
    } else if (any(familyName == MULTIVARIATE_FAMILY_NAMES)) {
      stop("Multivariate prior '", familyName, "' applied to coordinate '", coordinateName,
           "' at factor '", getFactorNameForNumber(regression, factorNumber), "'.");
    }
  }
  
  if (!any(ALL_FAMILY_NAMES == familyName)) {
    stop("Unrecognized or unsuitable parametric family name, \'", option, ".");
  }

  if (is.na(options)) options <- "";

  # evaluate the options string in the calling environment
  optionsList <- eval(parse(text=paste("list(", options, ")", sep="")),
                      envir=callingEnvironment);

  # Build up a list with all of the named options, filling in defaults where necessary.
  parsedPrior <- completeCorrelationPriorFromOptions(regression, callingEnvironment, factorNumber,
                                                     familyName, optionsList);

  if (any(UNIVARIATE_FAMILY_NAMES == familyName)) {
    if (isDefault) {
      result$coordinateDefault <- parsedPrior;
    } else {
      result$coordinatePriors[[coordinateNumber]] <- parsedPrior;
    }
  } else {
    result$matrixPrior <- parsedPrior;
  }
    
  return (result);
}

completeCorrelationPriorFromOptions <- function(regression, callingEnvironment, factorNumber, family, optionsList)
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
  fillDefaults <- function(x, ...) x;;

  if (family == GAMMA_FAMILY_NAME) {
    prior$shape <- getPriorOption(SHAPE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    rate  <- getPriorOption(RATE_HYPERPARAMETER_NAME,  namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    if (!is.null(rate) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$rate <- rate;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getCorrelationGammaDefaults;
  } else if (family == INVGAMMA_FAMILY_NAME) {
    prior$shape <- getPriorOption(SHAPE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    scale <- getPriorOption(SCALE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    if (!is.null(scale) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$scale <- scale;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getCorrelationInverseGammaDefaults;
  } else if (family == WISHART_FAMILY_NAME) {
    prior$degreesOfFreedom <-
      getPriorOption(DEGREES_OF_FREEDOM_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    scale <- getPriorOption(SCALE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    if (!is.null(scale) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$scale <- scale;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getCorrelationWishartDefaults;
  } else if (family == INVWISHART_FAMILY_NAME) {
    prior$degreesOfFreedom <-
      getPriorOption(DEGREES_OF_FREEDOM_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    inverseScale <- getPriorOption(INVERSE_SCALE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    if (!is.null(inverseScale) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$inverseScale <- inverseScale;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getCorrelationInverseWishartDefaults;
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


getCorrelationGammaDefaults <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (is.null(prior$shape))          prior$shape <- defaultCorrelationGammaShape;
  if (is.null(prior$posteriorScale)) prior$posteriorScale <- defaultCorrelationPosteriorScale;
  if (is.null(prior$rate))           prior$rate <- defaultCorrelationGammaRate;
  if (is.null(prior$dataScale))      prior$dataScale <- defaultCorrelationDataScale;

  return(evaluateGammaExpressions(regression, callingEnvironment, factorNumber, prior));
}

getCorrelationInverseGammaDefaults <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (is.null(prior$shape))          prior$shape <- defaultCorrelationInverseGammaShape;
  if (is.null(prior$posteriorScale)) prior$posteriorScale <- defaultCorrelationPosteriorScale;
  if (is.null(prior$scale))          prior$scale <- defaultCorrelationInverseGammaScale;
  if (is.null(prior$dataScale))      prior$dataScale <- defaultCorrelationDataScale;

  return(evaluateInverseGammaExpressions(regression, callingEnvironment, factorNumber, prior));
}

getCorrelationWishartDefaults <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (is.null(prior$degreesOfFreedom))
    prior$degreesOfFreedom <- defaultCorrelationWishartDegreesOfFreedom;
  if (is.null(prior$scale))
    prior$scale <- defaultCorrelationWishartScale;
  if (is.null(prior$dataScale))
    prior$dataScale <- defaultCorrelationDataScale;
  
  return(evaluateWishartExpressions(regression, callingEnvironment, factorNumber, prior));
}

getCorrelationInverseWishartDefaults <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (is.null(prior$degreesOfFreedom))
    prior$degreesOfFreedom <- defaultCorrelationInverseWishartDegreesOfFreedom;
  if (is.null(prior$inverseScale))
    prior$inverseScale <- defaultCorrelationInverseWishartInverseScale;
  if (is.null(prior$dataScale))
    prior$dataScale <- defaultCorrelationDataScale;
  
  return(evaluateInverseWishartExpressions(regression, callingEnvironment, factorNumber, prior));
}
