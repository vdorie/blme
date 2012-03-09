parseUnmodeledCoefficientPriorSpecification <-
  function(regression, specification, callingEnvironment)
{
  if (is.null(specification)) return(createFlatPriorObject());

  if (!is.character(specification)) {
    stop("Unmodeled coefficient prior specification is not of character type ",
         "and is unsupported.");
  }

  if (!grepl(unmodeledCoefficientPriorSpecificationPattern, specification, perl=TRUE)) {
    stop(paste("\"", specification, "\" does not parse as a prior specification", sep=""));
  }

  parse <- subSplit(unmodeledCoefficientPriorSpecificationPattern, "\\1\1\\2",
                    specification, perl=TRUE);
  priorType <- parse[1];
  options   <- parse[2];
  if (is.na(options)) options <- "";

  if (priorType == NONE_TYPE_NAME || priorType == FLAT_FAMILY_NAME)
    return (createFlatPriorObject());

  if (priorType == NORMAL_FAMILY_NAME) {
    prior <- parseUnmodeledCoefficientNormalPrior(regression, options, callingEnvironment);
  } else if (priorType == MVT_FAMILY_NAME) {
    prior <- parseUnmodeledCoefficientMVTPrior(regression, options, callingEnvironment);
  } else {
    stop("Internal error, please contact the package authors: prior type '",
         priorType, "' unsupported, yet parsed successfully.");
  }

  prior <- validateUnmodeledCoefficientPrior(regression, prior);
  prior <- adjustUnmodeledCoefficientPriorScales(regression, prior);

  fields <- getUnmodeledCoefficientPriorFields(regression, prior);
  
  return(createPriorObject(DIRECT_TYPE_NAME, fields));
}

parseUnmodeledCoefficientMVTPrior <- function(regression, options, callingEnvironment)
{
  stop("Multivariate-t distribution not yet implemented.");
}


getUnmodeledCoefficientDefault <- function(regression, family)
{
  result <- list(family = family);

  fillDefaults <- function(x, ...) x;
  if (result$family == NORMAL_FAMILY_NAME) {
    fillDefaults <- getUnmodeledCoefficientNormalDefaults;
  }

  result <- fillDefaults(result);
  return(result);
}

parseUnmodeledCoefficientNormalPrior <- function(regression, specification, callingEnvironment)
{
  errorPrefix <- "Error applying prior to unmodeled coefficients: ";
  numUnmodeledCoefs <- regression@dims[["p"]];
  
  option <- trim(specification);
  if (nchar(option) == 0) return(getUnmodeledCoefficientDefault(regression, NORMAL_FAMILY_NAME));
  
  optionsList <- eval(parse(text=paste("list(", option, ")", sep="")),
                      envir=callingEnvironment);
  
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

  family <- NORMAL_FAMILY_NAME; # would otherwise be passed in, but nothing else is supported
  
  fillDefaults <- function(x, ...) x;
  if (family == NORMAL_FAMILY_NAME) {
    prior <- list(family = NORMAL_FAMILY_NAME);

    standardDeviation <- getPriorOption(SD_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    covariance <- NULL;
    if (is.null(standardDeviation)) {
      covariance <- getPriorOption(COVARIANCE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    } else {
      covariance <- standardDeviation^2;
    }
    prior$covarianceScale <- getPriorOption(COVARIANCE_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
    dataScale <- getPriorOption(DATA_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);
  
    if (!is.null(covariance) && is.null(dataScale)) dataScale <- ABSOLUTE_SCALE_NAME;
    prior[["covariance"]] <- covariance;
    prior$dataScale <- dataScale;

    fillDefaults <- getUnmodeledCoefficientNormalDefaults;
  }

  if (length(namedOptions) > 0) {
    warning("Unrecognized prior option(s) for ", family, " family: ",
            toString(names(namedOptions)), ".");
  }
  if (length(unnamedOptions) > 0) {
    warning("Extra option(s) for ", family, " family: ", toString(unnamedOptions), ".");
  }
  
  prior <- fillDefaults(prior);

  return(prior);
}

isDiagonal <- function(matrix) {
  if (!is.matrix(matrix)) return(FALSE);
  if (nrow(matrix) != ncol(matrix)) return(FALSE);
  if (nrow(matrix) == 1) return(TRUE);
  
  return(all(matrix[upper.tri(matrix)] == 0) &&
         all(matrix[lower.tri(matrix)] == 0));
}

getUnmodeledCoefficientPriorFields <- function(regression, prior)
{
  numUnmodeledCoefs <- regression@dims[["p"]];
  
  families <- integer(0);
  scales   <- integer(0);
  hyperparameters <- double(0);
  
  if (prior$family == NORMAL_FAMILY_NAME) {
    families <- getEnumOrder(familyEnumeration, NORMAL_FAMILY_NAME);
    scales   <- getEnumOrder(scaleEnumeration, prior[["covarianceScale"]]);

    if (length(scales) == 0) {
      stop("Unable to recognize covariance scale '", prior[["covarianceScale"]], "'.");
    }

    covariance <- prior[["covariance"]];
    if (length(covariance) == 1) {
      logCovDet <- numUnmodeledCoefs * log(covariance);
      hyperparameters <- c(logCovDet, 1 / sqrt(covariance));
    } else if (length(covariance) == numUnmodeledCoefs) {
      logCovDet <- sum(log(covariance));
      hyperparameters <- c(logCovDet, 1 / sqrt(covariance));
    } else if (isDiagonal(covariance)) {
      logCovDet <- sum(log(diag(covariance)));
      hyperparameters <- c(logCovDet, 1 / sqrt(diag(covariance)));
    } else {
      covariance.inv <- solve(covariance);
      leftFactor <- solve(t(chol(covariance)));
      covarianceLength <- numUnmodeledCoefs^2;

      logCovDet <- determinant(covariance, logarithm=TRUE)$modulus;
      hyperparameters <- rep(0, 1 + 2 * covarianceLength)
      hyperparameters[1] <- logCovDet;
      hyperparameters[2:(covarianceLength + 1)] <- leftFactor;
      hyperparameters[(covarianceLength + 2):(2 * covarianceLength + 1)] <- covariance.inv;
    }
  }

  return (list(families        = families,
               scales          = scales,
               hyperparameters = hyperparameters));
}

getUnmodeledCoefficientNormalDefaults <- function(prior)
{
  if (is.null(prior[["covarianceScale"]]))
    prior[["covarianceScale"]] <- defaultUnmodeledCoefficientNormalCovarianceScale;
  if (is.null(prior[["covariance"]]))
    prior[["covariance"]] <- defaultUnmodeledCoefficientNormalSD^2;
  if (is.null(prior$dataScale))
    prior$dataScale <- defaultUnmodeledCoefficientDataScale;

  return(prior);
}

unmodeledCoefficientPriorToString <- function(regression)
{
  if (regression@fixef.prior@type == getEnumOrder(typeEnumeration, NONE_TYPE_NAME))
    return(character(0));

  
  families <- regression@fixef.prior@families;
  scales <- regression@fixef.prior@scales;
  hyperparameters <- regression@fixef.prior@hyperparameters;

  numUnmodeledCoefficients <- regression@dims[["p"]];
  if (length(hyperparameters) == 1 + 1) {
    # just a single sd
    hyperparameters <- 1 / hyperparameters[2]^2;
  } else if (length(hyperparameters) == 1 + numUnmodeledCoefficients) {
    # diagonal of sd
    hyperparameters <- 1 / hyperparameters[1 + 1:numUnmodeledCoefficients]^2;
  } else {
    # full matrix
    covarianceInverse <- matrix(hyperparameters[1 + numUnmodeledCoefficients^2 +
                                                1:(numUnmodeledCoefficients^2)],
                                numUnmodeledCoefficients,
                                numUnmodeledCoefficients);
    hyperparameters <- as.numeric(solve(covarianceInverse));
  }

  return(buildStringForFamily(families, scales, hyperparameters, TRUE)$string);
}
