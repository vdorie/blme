getSpectralDefault <- function(regression, callingEnvironment, factorNumber)
{
  result <- list(family = defaultSpectralFamily);

  fillDefaults <- function(x, ...) x;
  if (result$family == GAMMA_FAMILY_NAME) {
    fillDefaults <- getSpectralGammaDefaults;
  } else if (result$family == INVGAMMA_FAMILY_NAME) {
    fillDefaults <- getSpectralInverseGammaDefaults;
  }
  
  prior <- fillDefaults(regression, callingEnvironment, factorNumber, result);
  if (factorNumber > 0) return(prior);
  
  return(list(prior = prior, fillDefaults = fillDefaults));
}

parseSpectralSpecification <- function(regression, callingEnvironment, factorNumber,
                                       specification)
{  
  if (factorNumber == 0) {
    priors <- list();
  } else {
    factorDimension <- nrow(regression@ST[[factorNumber]]);
    priors          <- vector("list", factorDimension);
  }
  
  default <- getSpectralDefault(regression, callingEnvironment, factorNumber);
  
  option <- trim(specification);
  while (nchar(option) > 0) {
    # should go through the options one at a time, identifying any defaults
    # and installing the right priors in the right spots
    parseResult <- parseSpectralOption(regression, callingEnvironment, factorNumber,
                                       priors, default, option);
    
    priors  <- parseResult$priors;
    default <- parseResult$default;
    option  <- parseResult$option;
  }

  # convert the result to something that the C can use sensibly
  result <- NULL;
  if (factorNumber == 0) {
    result <- default;
    result$type <- SPECTRAL_TYPE_NAME;
  } else {
    result <- list(type = SPECTRAL_TYPE_NAME, prior = priors, coordinateDefault = default);
  }

  return (result);
}


parseSpectralOption <- function(regression, callingEnvironment, factorNumber,
                                priors, default, option)
{
  # split the first coordinate prior off from the rest of the comma separated list
  commaSplit <- splitListOnComma(option);
  if (length(commaSplit) == 1) {
    remainder <- "";
  } else {
    option    <- commaSplit[1];
    remainder <- trim(commaSplit[2]);
  }

  result <- list(priors = priors, default = default, option = remainder);

  supportedFamilies <- flattenStrings(UNIVARIATE_FAMILY_NAMES, sep="|"); # for use in regular expressions
  
  isDefault <- FALSE;
  coordinateNumber <- 0;
  # figure out the coordinate name (or none if default), parametric family, and the options
  #
  # if the factorNumber is 0, i.e. is a global default, for now we ignore any coordinate specifics
  if (isDefaultSpecification(option) || factorNumber == 0) {
    splitFamilyFromOptionsRegex <- paste("^(", supportedFamilies, ")",
                                         "(?:\\s*\\((.*)\\))?$", sep="");
    
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

    coordinateNumber <- tryCatch(as.integer(coordinateName), warning = function(e) NA);
    if (is.na(coordinateNumber)) {
      stop("Unable to identify coordinate number in '", option, "'.");
    }
    if (coordinateNumber < 1 || coordinateNumber > nrow(regression@ST[[factorNumber]]))
      stop("Coordinate number ", coordinateNumber, " out of range.");
  }
  
  if (!any(UNIVARIATE_FAMILY_NAMES == familyName)) {
    stop("Unrecognized or unsuitable parametric family name, '", option, "'.");
  }

  if (is.na(options)) options <- "";

  # evaluate the options string in the calling environment
  optionsList <- eval(parse(text=paste("list(", options, ")", sep="")),
                      envir=callingEnvironment);

  parsedPrior <- completeSpectralPriorFromOptions(regression, callingEnvironment, factorNumber,
                                                  familyName, optionsList);
 
  if (isDefault) {
    result$default <- parsedPrior;
  } else {
    result$priors[[coordinateNumber]] <- parsedPrior;
  }
  
  return (result);
}


completeSpectralPriorFromOptions <- function(regression, callingEnvironment, factorNumber,
                                             family, optionsList)
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
    rate  <- getPriorOption(RATE_HYPERPARAMETER_NAME,  namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    if (!is.null(rate) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$rate <- rate;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getDirectGammaDefaults;
  } else if (family == INVGAMMA_FAMILY_NAME) {
    prior$shape <- getPriorOption(SHAPE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    scale <- getPriorOption(SCALE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    if (!is.null(scale) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior$scale <- scale;
    prior$dataScale <- dataScale;
    
    fillDefaults <- getDirectInverseGammaDefaults;
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


getSpectralGammaDefaults <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (is.null(prior$shape))          prior$shape <- defaultSpectralGammaShape;
  if (is.null(prior$posteriorScale)) prior$posteriorScale <- defaultSpectralPosteriorScale;
  if (is.null(prior$rate))           prior$rate <- defaultSpectralGammaRate;
  if (is.null(prior$dataScale))      prior$dataScale <- defaultSpectralDataScale;

  return(evaluateGammaExpressions(regression, callingEnvironment, factorNumber, prior));
}

getSpectralInverseGammaDefaults <- function(regression, callingEnvironment, factorNumber, prior)
{
  if (is.null(prior$shape))          prior$shape <- defaultSpectralInverseGammaShape;
  if (is.null(prior$posteriorScale)) prior$posteriorScale <- defaultSpectralPosteriorScale;
  if (is.null(prior$scale))          prior$scale <- defaultSpectralInverseGammaScale;
  if (is.null(prior$dataScale))      prior$dataScale <- defaultSpectralDataScale;

  return(evaluateInverseGammaExpressions(regression, callingEnvironment, factorNumber, prior));
}
