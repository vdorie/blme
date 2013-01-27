# the prior specification may look like:
#   a) a single character string
#   b) a list of strings of the form in a
#
# any single string looks like:
#
#   factorName ~ priorType(coordinatePrior1, coordinatePrior2, ..., [matrixPrior])
#
#
# factorName should match one of the dimnames of the factor list
# if factorName is omitted, a default prior is set
#
# priorType can be "flat", "spectral", "correlation", or any of a
# supported set of distributions
#
#   coordinate priors have the form:
#
#     coordinateName ~ family(posterior.scale = scale, ...)
#
#
#   the coordinate name can be omitted to set a default. numbers starting at 1
#   are taken to specify the coordinate by index (as only makes sense in spectral
#   case).
#
#   posterior.scale should be "sd" or "var". Default is "sd" for the correlation prior and
#   "var" for the spectral prior. In the case of spectral priors, "sd" scale is the square
#   root of the eigenvalues.
#
#   ... should be a list of parameters
#
#   supported families include:
#     for spectral and correlation
#       gamma(shape=, rate=)
#       inverse.gamma(shape=, scale=)
#     directly applied
#       gamma(shape=, rate=, posterior.scale) (on univariate factors)
#       inverse.gamma(shape=, scale=, posterior.scale) (on univariate factors)
#       wishart(df=,scale=)
#       inverse.wishart(df=,inverse.scale=)
#
#
#   the matrix prior currently only applies to the correlation prior type, and
#   has form:
#     family(...)
#
#   supported families include:
#     wishart(df = degreesOfFreedom, scale = matrix)
#     inverse.wishart(df=,inverse.scale=)


parseCovariancePriorSpecification <- function(regression, priorSpecification, callingEnvironment)
{
  if (is.null(priorSpecification)) return (getDefaultCovariancePrior(regression));
  
  if (!is.list(priorSpecification) && !is.character(priorSpecification)) {
    stop("Covariance prior specification is not of character or list type ",
         "and is unsupported.");
  }

  if (!is.list(priorSpecification)) {
    # turns a comma separated list into an actual list, respecting parentheses
    priorSpecification <- trim(splitListOnCommas(priorSpecification));
  }
  
  # Keep a running copy of the default, and start out with a "NULL" prior.
  # NULL prior differs from default in that it is all set to "nothingness".
  # If any of the nulls are replaced, we know not to clobber them with
  # a default.
  numFactors <- regression@dims[["nt"]];
  
  default <- list(univariate = list(), multivariate = list());
  prior   <- vector("list", numFactors);

  for (i in 1:length(priorSpecification)) {
    # parseSpecification pops off the head of the list and installs it
    # as a default or at a particular grouping factor.
    parseResult <- parseCovarianceSpecification(regression, callingEnvironment,
                                                priorSpecification[[i]],
                                                prior, default);
    
    prior   <- parseResult$prior;
    default <- parseResult$default;
  }
  
  # install defaults if necessary
  for (factorNumber in 1:numFactors) {
    factorPrior <- prior[[factorNumber]];
    if (is.null(factorPrior)) {
      factorPrior <- createDefaultPriorForFactor(regression, callingEnvironment, factorNumber, default);
    }
    
    if (factorPrior$type != NONE_TYPE_NAME) {
      factorPrior$prior <- validateCovariancePrior(regression, factorPrior$prior, factorNumber);
      factorPrior$prior <- adjustCovariancePriorScales(regression, factorPrior$prior, factorPrior$type,
                                                     factorNumber);
    
      factorDimension <- nrow(regression@ST[[factorNumber]]);
      fields <- getPriorFields(factorDimension, factorPrior$prior, factorPrior$coordinateDefault);
      prior[[factorNumber]] <- createPriorObject(factorPrior$type, fields);
    } else {
      prior[[factorNumber]] <- createFlatPriorObject();
    }
  }
  
  return (prior);
}

createDefaultPriorForFactor <- function(regression, callingEnvironment, factorNumber, default) {
  factorDimension <- nrow(regression@ST[[factorNumber]]);
  
  # For univariate factors, we might have a special default.
  if (factorDimension == 1) {
    if (!is.null(default$univariate) && length(default$univariate) > 0 &&
        default$univariate$type != NONE_TYPE_NAME) {
      prior <- getDefaultUnivariatePriorForFactor(regression, callingEnvironment,
                                                  factorNumber, default$univariate);
      return(list(type = DIRECT_TYPE_NAME,
                  prior = prior));
    }
    
    if (is.null(default$multivariate) || length(default$multivariate) == 0)
      return(list(type = NONE_TYPE_NAME));
  }
  
  # If we're still going, we can assume that the mv prior should apply.
  # Throws out (for this function call) the univariate default.
  default <- default$multivariate;
  if (length(default) == 0) return(list(type = NONE_TYPE_NAME));
  
  if (default$type == SPECTRAL_TYPE_NAME) {
    # a spectral prior is a univariate one, applied multiple times
    coordinatePriors <- getDefaultUnivariatePriorForFactor(regression, callingEnvironment,
                                                           factorNumber, default);
    return(list(type = SPECTRAL_TYPE_NAME,
                prior = coordinatePriors));
    
  } else if (default$type == CORRELATION_TYPE_NAME) {
    univariatePriors <- getDefaultUnivariatePriorForFactor(regression, callingEnvironment,
                                                           factorNumber, default$coordinate);
    multivariatePrior <- getDefaultMultivariatePriorForFactor(regression, callingEnvironment,
                                                              factorNumber, default$matrix);
    
    return(list(type = CORRELATION_TYPE_NAME,
                prior = append(univariatePriors, list(multivariatePrior))));
  } else if (default$type == DIRECT_TYPE_NAME) {
    multivariatePrior <- list(getDefaultMultivariatePriorForFactor(regression, callingEnvironment,
                                                                   factorNumber, default));
    return(list(type = DIRECT_TYPE_NAME,
                prior = multivariatePrior));
  } else if (default$type == NONE_TYPE_NAME) {
    return(list(type = NONE_TYPE_NAME));
  }
  
  # shouldn't happen, as an error should have been raised
  # when we were unable to parse this field
  stop("Prior of type '", default$type, "' unrecognized.");

  return (createFlatPriorObject());
}

getDefaultUnivariatePriorForFactor <- function(regression, callingEnvironment,
                                               factorNumber, default)
{
  prior <- default$prior;
  fillDefaults <- default$fillDefaults;

  result <- prior;
  result <- fillDefaults(regression, callingEnvironment, factorNumber, result);

  factorDimension <- nrow(regression@ST[[factorNumber]]);
  
  return(lapply(1:factorDimension, function(i) result));
}

getDefaultMultivariatePriorForFactor <- function(regression, callingEnvironment,
                                                 factorNumber, default)
{
  prior <- default$prior;
  fillDefaults <- default$fillDefaults;

  result <- prior;
  result <- fillDefaults(regression, callingEnvironment, factorNumber, result);

  return(result);
}

mergePriorFields <- function(fields1, fields2) {
  allFieldNames <- union(names(fields1), names(fields2));
  result <- list();
  for (fieldName in allFieldNames) {
    if (is.null(fields1[[fieldName]])) {
      result[[fieldName]] <- fields2[[fieldName]];
    } else if (is.null(fields2[[fieldName]])) {
      result[[fieldName]] <- fields1[[fieldName]];
    } else {
      result[[fieldName]] <- append(fields1[[fieldName]], fields2[[fieldName]]);
    }
  }
      
  return (result);
}

# installs a flat prior on everything 
getDefaultCovariancePrior <- function(regression)
{
  numFactors = regression@dims[["nt"]];
  prior <- vector("list", numFactors);
  
  for (i in 1:numFactors) {
    prior[[i]] <- list(type = NONE_TYPE_NAME);
  }
  return (prior);
}

# When this function is called, we have something like:
#
#   factorName ~ priorType
#   factorName ~ priorType(optionsList)
#   defaultPriorType
#   defaultPriorType(optionsList)
parseCovarianceSpecification <- function(regression, callingEnvironment, specification, prior, default) {
  result <- list(prior = prior, default = default);
  
  if (is.null(specification) || !is.character(specification)) {
    stop("Component of covariance prior not of character type.");
  }

  isValidSpecification <-
    grepl(covariancePriorSpecificationPattern, specification, perl=TRUE);
  if (!isValidSpecification) {
    stop("'", specification, "' does not parse as a prior specification.");
  }

  parse <- subSplit(covariancePriorSpecificationPattern, "\\1\1\\2\1\\3",
                    specification, perl=TRUE);
  
  factorName <- parse[1];
  priorType  <- parse[2];
  remainder  <- parse[3];
  if (is.na(remainder)) remainder <- "";

  isDefault <- FALSE;
  factorNumber <- tryCatch(as.integer(factorName), warning = function(e) NA);
  if (is.na(factorNumber)) {
    factorNumber <- getFactorNumberForName(regression, factorName);
    if (length(factorNumber) > 1) {
      warning("Factor name '", factorName, "' associated with multiple ",
              "terms. Only the first will be used.");
      factorNumber <- factorNumber[1];
    }
    if (factorNumber == 0 && factorName != "") {
      stop("Factor name '", factorName, "' unrecognized.");
    }
  }

  if (factorNumber < 0 || factorNumber > regression@dims[["nt"]]) {
    stop("Factor number ", factorNumber, " out of range.");
  }

  if (factorNumber == 0) isDefault <- TRUE;

  # supported 'types' are the "none", "flat", the decompositions, or any distributional family
  # should have failed before here and this test should be superfluous
  if (!any(priorType == SUPPORTED_COVARIANCE_TYPE_NAMES)) {
    stop("Unrecognized prior type, '", priorType, "'.");
  }

  if (priorType == NONE_TYPE_NAME || priorType == FLAT_FAMILY_NAME) {
    parse <- parseFlatSpecification(regression, factorNumber, remainder);

    # none could be univariate or multivariate
    if (isDefault && !is.null(parse)) {
      result$default$univariate   <- parse;
      result$default$multivariate <- parse;
    } else {
      result$prior[[factorNumber]] <- parse;
    }
  } else if (any(priorType == ALL_FAMILY_NAMES)) {
    # function defined in directPriorParser.R
    parse <- parseDirectSpecification(regression, callingEnvironment, factorNumber, priorType, remainder);
    
    if (isDefault && !is.null(parse)) {
      if (any(priorType == UNIVARIATE_FAMILY_NAMES)) {
        result$default$univariate <- parse;
      } else {
        result$default$multivariate <- parse;
      }
    } else {
      result$prior[[factorNumber]] <- parse;
    }
  } else {
    parse <- NULL;
    if (priorType == SPECTRAL_TYPE_NAME) {
      # function defined in spectralPriorParser.R
      parse <- parseSpectralSpecification(regression, callingEnvironment, factorNumber, remainder);
    } else if (priorType == CORRELATION_TYPE_NAME) {
      # function defined in correlationPriorParser.R
      parse <- parseCorrelationSpecification(regression, callingEnvironment, factorNumber, remainder);
    }

    if (isDefault && !is.null(parse)) {
      result$default$multivariate <- parse;
    } else {
      result$prior[[factorNumber]] <- parse;
    }
  }
  
  return(result);
}

parseFlatSpecification <- function(regression, factorNumber, specification) {
  return(list(type = NONE_TYPE_NAME));
}


# takes the prior in list form and returns it as vectors that the C side of things
# can easily use
getPriorFields <- function(factorDimension, priors, default)
{
  families <- integer(length(priors));
  
  numScales <- length(priors);
  numHyperparameters <- 0;
  matrixIndex   <- 0;   # as the matrix prior might not be last in the list, keep track

  # first count the length of vectors needed to create the prior object
  for (i in 1:length(priors)) {
    if (is.null(priors[[i]])) priors[[i]] <- default;
    
    if (priors[[i]]$family == GAMMA_FAMILY_NAME ||
        priors[[i]]$family == INVGAMMA_FAMILY_NAME)
    {
      numHyperparameters <- numHyperparameters + 2;
    } else if (priors[[i]]$family == WISHART_FAMILY_NAME ||
               priors[[i]]$family == INVWISHART_FAMILY_NAME)
    {
      # df + log(determinant) + entire scale matrix
      numHyperparameters <- numHyperparameters + 2 + factorDimension * factorDimension;
      matrixIndex <- i;
    }
  }
  scales <- integer(numScales);
  hyperparameters <- rep(0, numHyperparameters);

  # fill in the hyperparameters
  familyIndex <- 1;
  hyperparameterIndex <- 1;
  for (i in 1:length(priors)) {
    if (i == matrixIndex) next;
    
    families[familyIndex] <- getEnumOrder(familyEnum, priors[[i]]$family);
    posteriorScale <- getEnumOrder(posteriorScaleEnum, priors[[i]]$posteriorScale);
    commonScale    <- getEnumOrder(commonScaleEnum, priors[[i]]$commonScale);
    scales[familyIndex] <- getScaleInt(posteriorScale, commonScale);
    if (priors[[i]]$family == GAMMA_FAMILY_NAME) {
      hyperparameters[hyperparameterIndex]     <- priors[[i]]$shape;
      hyperparameters[hyperparameterIndex + 1] <- priors[[i]]$rate;
      hyperparameterIndex <- hyperparameterIndex + 2;
    } else if (priors[[i]]$family == INVGAMMA_FAMILY_NAME) {
      hyperparameters[hyperparameterIndex]     <- priors[[i]]$shape;
      hyperparameters[hyperparameterIndex + 1] <- priors[[i]]$scale;
      hyperparameterIndex <- hyperparameterIndex + 2;
    }
    familyIndex <- familyIndex + 1;
  }

  if (matrixIndex > 0) {
    posteriorScale <- getEnumOrder(posteriorScaleEnum, priors[[matrixIndex]]$posteriorScale);
    commonScale    <- getEnumOrder(commonScaleEnum, priors[[matrixIndex]]$commonScale);
    scales[matrixIndex] <- getScaleInt(posteriorScale, commonScale);
    # Order of wish parameters:
    #   df, log(det(scale)), scale
    # iwish is similar, but some inverses are required
    if (priors[[matrixIndex]]$family == INVWISHART_FAMILY_NAME) {
      inverseScaleMatrix <- priors[[matrixIndex]]$inverseScale;
      logInverseScaleDet <- determinant(inverseScaleMatrix, logarithm=TRUE)$modulus;

      hyperparameters[hyperparameterIndex] <- priors[[matrixIndex]]$degreesOfFreedom;
      hyperparameters[hyperparameterIndex + 1] <- logInverseScaleDet;
      hyperparameterIndex <- hyperparameterIndex + 2;
      
      hyperparameters[hyperparameterIndex:(hyperparameterIndex + factorDimension * factorDimension - 1)] <-
        inverseScaleMatrix;
      hyperparameterIndex <- hyperparameterIndex + factorDimension * factorDimension;      
    } else if (priors[[matrixIndex]]$family == WISHART_FAMILY_NAME) {
      scaleMatrix <- priors[[matrixIndex]]$scale;
      logScaleDet <- determinant(scaleMatrix, logarithm=TRUE)$modulus;
      if (is.finite(logScaleDet)) {
        scaleInverse <- as.vector(solve(scaleMatrix));
      } else {
        scaleInverse <- rep(0, length(scaleMatrix));
      }
    
      hyperparameters[hyperparameterIndex] <- priors[[matrixIndex]]$degreesOfFreedom;
      hyperparameters[hyperparameterIndex + 1] <- logScaleDet;
      hyperparameterIndex <- hyperparameterIndex + 2;
      
      hyperparameters[hyperparameterIndex:(hyperparameterIndex + factorDimension * factorDimension - 1)] <-
        scaleInverse;
      hyperparameterIndex <- hyperparameterIndex + factorDimension * factorDimension;
    }
      
    families[familyIndex] <- getEnumOrder(familyEnum, priors[[matrixIndex]]$family);
    familyIndex <- familyIndex + 1;
  }

  return (list(families        = families,
               scales          = scales,
               hyperparameters = hyperparameters));
}

covariancePriorToString <- function(regression)
{
  stringResult <- NULL; # for R CMD Check
  
  stringConnection <- textConnection("stringResult", "w", local=TRUE);
  sink(stringConnection);
  
  numFactors <- regression@dims[["nt"]];
  factorNames <- expandFactorNames(regression);
  for (i in 1:numFactors) {
    prior <- regression@cov.prior[[i]];
    
    families <- prior@families;
    scales <- prior@scales;
    hyperparameters <- prior@hyperparameters;
    
    if (prior@type == getEnumOrder(typeEnum, NONE_TYPE_NAME)) next;
    
    cat(factorNames[i], "~ ");
    
    if (prior@type == getEnumOrder(typeEnum, DIRECT_TYPE_NAME)) {
      cat(buildStringForFamily(families, scales, hyperparameters, 2, TRUE)$string,
          "\n", sep = "");
    } else if (prior@type == getEnumOrder(typeEnum, CORRELATION_TYPE_NAME)) {
      coordinateNames <- colnames(regression@ST[[i]]);
      
      cat(typeEnum[prior@type + 1], "\n", sep = "");
      
      for (j in 1:length(coordinateNames)) {
        familyString <- buildStringForFamily(families, scales, hyperparameters, 2, TRUE);
        cat("  ", coordinateNames[j], " ~ ", stringResult$familyString, "\n", sep = "");

        families <- families[(familyString$numFamilies + 1):length(families)];
        scales   <- scales[(familyString$numScales + 1):length(families)];
        hyperparameters <- hyperparameters[(familyString$numHyperparameters + 1):
                                           length(hyperparameters)];
        
      }
      cat("  ", buildStringForFamily(families, scales, hyperparameters, 2, TRUE), "\n", sep="");
    } else if (prior@type == getEnumOrder(typeEnum, SPECTRAL_TYPE_NAME)) {
      coordinateNames <- colnames(regression@ST[[i]]);
      
      cat(typeEnum[prior@type + 1], "\n", sep = "");
      
      for (j in 1:length(coordinateNames)) {
        familyString <- buildStringForFamily(families, scales, hyperparameters, 2, TRUE);
        cat("  ", j, " ~ ", familyString$string, "\n", sep = "");
        
        families <- families[(familyString$numFamilies + 1):length(families)];
        scales   <- scales[(familyString$numScales + 1):length(families)];
        hyperparameters <- hyperparameters[(familyString$numHyperparameters + 1):
                                           length(hyperparameters)];
      }
    } else cat("unknown\n");
  }
  
  sink();
  close(stringConnection);

  return(stringResult);
}
