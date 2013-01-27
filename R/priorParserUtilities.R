# functions that the prior parser might call that don't influence the logic of the
# of the parse
getEnumOrder <- function(Enum, name)
{
  return (as.integer(which(Enum == name) - 1L));
}

getScaleInt <- function(posteriorScale, commonScale) {
  if (is.na(posteriorScale) || is.na(commonScale) ||
      !is.numeric(posteriorScale) || !is.numeric(commonScale) ||
      (posteriorScale != 0L && posteriorScale != 1L) ||
      (commonScale != 0L && commonScale != 1L)) {
    stop("blme internal error: getScaleInt called on non-binary inputs ('",
         posteriorScale, "' and '", commonScale, "').");
  }
  return(.Call(bmer_getScaleInt, as.integer(posteriorScale), as.integer(commonScale)));
}

isDefaultSpecification <- function(specification) {
  if (is.null(specification) || typeof(specification) != "character") return (FALSE);

  patternMatches <- grepl(defaultSpecificationRegularExpression, specification, perl=TRUE)

  return (patternMatches);
}

expandFactorNames <- function(regression)
{
  names(regression@flist)[attr(regression@flist, "assign")];
}

getFactorNumberForName <- function(regression, factorName)
{
  factorNames <- expandFactorNames(regression);
  if (!any(factorNames == factorName)) return(0);
  return (which(factorNames == factorName));
}

getFactorNameForNumber <- function(regression, factorNumber)
{
  return(expandFactorNames(regression)[factorNumber]);
}

getCoordinateNumberForName <- function(regression, factorNumber, coordinateName)
{
  coordinateNames <- colnames(regression@ST[[factorNumber]]);
  if (!any(coordinateNames == coordinateName)) return(0);
  return (which(coordinateNames == coordinateName));
}

getInputSdForUnmodeledCoefficients <- function(regression)
{
  numUnmodeledCoef <- regression@dims[["p"]];
  result <- rep(1, numUnmodeledCoef);

  if (numUnmodeledCoef == 1) return(result);
  
  for (i in 2:numUnmodeledCoef) {
    result[i] <- sd(regression@X[,i]);
  }
    
#  columnNames <- names(regression@fixef);
#  result <- rep(NA, numUnmodeledCoef);
#  for (i in 1:numUnmodeledCoef) {
#    if (columnNames[i] == "(Intercept)") {
#      result[i] <- 1;
#    } else {
#      result[i] <- sd(regression@frame[[columnNames[i]]]);
#      if (result[i] == 0) result[i] <- 1;
#    }
#  }

  return(result);
}

getInputSdForFactor <- function(regression, factorNumber)
{
  factorDimension <- nrow(regression@ST[[factorNumber]]);
  factorPredictorNames <- rownames(regression@ST[[factorNumber]]);
  fixefNames <- names(regression@fixef);    # fixef matches with X matrix
  frameNames <- colnames(regression@frame); # in case it's in the frame but not a fixef 

  result <- rep(NA, factorDimension);
  for (i in 1:factorDimension) {
    predictorName <- factorPredictorNames[i];

    if (predictorName == "(Intercept)") {
      result[i] <- 1;
    } else {
      if (any(fixefNames == predictorName)) {
        colNum <- which(fixefNames == predictorName);
        result[i] <- sd(regression@X[,colNum]);
      } else if (any(frameNames == predictorName)) {
        colNum <- which(frameNames == predictorName);
        result[i] <- sd(regression@frame[,colNum]);
      } else {
        warning("Unable to standardize input '", predictorName,
                "' for factor number ", factorNumber, ".");
      }
      if (is.na(result[i]) || result[i] == 0) result[i] <- 1;
    }
  }

  return(result);
}

# Function behaves as follows:
#   If the parameter is can be found in the named list
#     Return that and remove named element from list
#   If not and there are unnamed options
#     Return the first item in unnamedOptions and remove it from list
#   Otherwise, return default
#
# Achieves the removal of elements from the parent from by having passed in
# "references", which consist of the environment for which the named
# variable should be updated.
getPriorOption <- function(name, namedOptionsReference, unnamedOptionsReference)
{
  namedOptionsEnv <- namedOptionsReference$env;
  namedOptions    <- namedOptionsReference$var;
  unnamedOptionsEnv <- unnamedOptionsReference$env;
  unnamedOptions    <- unnamedOptionsReference$var;

  result <- NULL;
  if (!is.null(namedOptionsEnv[[namedOptions]][[name]])) {
    result <- namedOptionsEnv[[namedOptions]][[name]];

    currentList <- namedOptionsEnv[[namedOptions]];
    namedOptionsEnv[[namedOptions]] <- currentList[names(currentList) != name];
  } else if (length(unnamedOptionsEnv[[unnamedOptions]]) >= 1) {
    result <- unnamedOptionsEnv[[unnamedOptions]][[1]];

    currentList <- unnamedOptionsEnv[[unnamedOptions]];
    if (length(currentList) == 1) {
      unnamedOptionsEnv[[unnamedOptions]] <- list();
    } else {
      unnamedOptionsEnv[[unnamedOptions]] <- currentList[2:length(currentList)]
    }
  }

  return (result);
}


isLinearMixedModel <- function(regression) {
  return (length(regression@muEta) == 0 &&
          length(regression@V) == 0);
}


buildStringForFamily <- function(families, scales, hyperparameters, digits, preprocessed)
{
  matchedCall <- match.call();
  
  numFamiliesUsed <- 0;
  numScalesUsed <- 0;
  numHyperparametersUsed <- 0;

  stringResult <- NULL; # for R CMD Check
  
  stringConnection <- textConnection("stringResult", "w", local=TRUE);
  sink(stringConnection);

  cat(familyEnum[families[1] + 1]);

  # TODO: make this not hard-coded; can't seem to find any bit-wise & or |
  posteriorScale <- ifelse(scales[1] == 1 || scales[1] == 3, 1, 0) + 1;
  commonScale    <- ifelse(scales[1] == 2 || scales[1] == 3, 1, 0) + 1;

  if (families[1] == getEnumOrder(familyEnum, GAMMA_FAMILY_NAME)) {
    hyperparameters <- round(hyperparameters, digits);
    cat("(", SHAPE_HYPERPARAMETER_NAME, " = ", hyperparameters[1],
        ", ", RATE_HYPERPARAMETER_NAME, " = ", hyperparameters[2],
        ", ", POSTERIOR_SCALE_OPTION_NAME, " = ", posteriorScaleEnum[posteriorScale],
        ", ", COMMON_SCALE_OPTION_NAME, " = ", commonScaleEnum[commonScale],
        ")", sep="");
    numFamiliesUsed <- 1;
    numScalesUsed <- 1;
    numHyperparametersUsed <- 2;
  } else if (families[1] == getEnumOrder(familyEnum, INVGAMMA_FAMILY_NAME)) {
    hyperparameters <- round(hyperparameters, digits);
    cat("(", SHAPE_HYPERPARAMETER_NAME, " = ", hyperparameters[1],
        ", ", SCALE_HYPERPARAMETER_NAME, " = ", hyperparameters[2],
        ", ", POSTERIOR_SCALE_OPTION_NAME, " = ", posteriorScaleEnum[posteriorScale],
        ", ", COMMON_SCALE_OPTION_NAME, " = ", commonScaleEnum[commonScale],
        ")", sep="");
    numFamiliesUsed <- 1;
    numScalesUsed <- 1;
    numHyperparametersUsed <- 2;
  } else if (families[1] == getEnumOrder(familyEnum, WISHART_FAMILY_NAME)) {
    logScaleDet <- hyperparameters[2];
    scaleInverse <- hyperparameters[3:length(hyperparameters)];
    hyperparameters <- round(hyperparameters, digits);
    levelDimension <- round(sqrt(length(scaleInverse)), digits=0);
    scaleInverse <- matrix(scaleInverse, levelDimension, levelDimension);
    
    if (is.finite(logScaleDet)) {
      scale <- as.numeric(solve(scaleInverse));
    } else {
      scale <- matrix(Inf, levelDimension, levelDimension);
    }
    
    cat("(", DEGREES_OF_FREEDOM_HYPERPARAMETER_NAME, " = ", hyperparameters[1], ", ",
        SCALE_HYPERPARAMETER_NAME, " = ", sep="");
    if (levelDimension == 1) {
      cat(round(scale[1], digits));
    } else {
      if (levelDimension > 2)
        cat("c(", toString(format(scale[1:4], digits=digits, scientific=TRUE)), ", ...)", sep="")
      else
        cat("c(", toString(format(scale, digits=digits, scientific=TRUE)), ")", sep="");
    }
    cat(", ", COMMON_SCALE_OPTION_NAME, " = ", commonScaleEnum[commonScale], sep="");
    cat(")");
  } else if (families[1] == getEnumOrder(familyEnum, INVWISHART_FAMILY_NAME)) {
    hyperparameters <- round(hyperparameters, digits);
    inverseScale <- hyperparameters[3:length(hyperparameters)];
    levelDimension <- round(sqrt(length(inverseScale)), digits=0);
    
    cat("(", DEGREES_OF_FREEDOM_HYPERPARAMETER_NAME, " = ", hyperparameters[1], ", ",
        INVERSE_SCALE_HYPERPARAMETER_NAME, " = ", sep="");
    if (levelDimension == 1) {
      cat(inverseScale[1]);
    } else {
      if (levelDimension > 2)
        cat("c(", toString(format(inverseScale[1:4], digits=digits, scientific=TRUE)), ", ...)", sep="")
      else
        cat("c(", toString(format(inverseScale, digits=digits, scientific=TRUE)), ")", sep="");
    }
    cat(", ", COMMON_SCALE_OPTION_NAME, " = ", commonScaleEnum[commonScale], sep="");
    cat(")");
  } else if (families[1] == getEnumOrder(familyEnum, NORMAL_FAMILY_NAME)) {
    hyperparameters <- round(hyperparameters, digits);
    numParams <- length(hyperparameters);
    
    if (preprocessed) {
      cat("(", COVARIANCE_HYPERPARAMETER_NAME, " = ", sep = "");
    } else {
      cat("(hyperparams = ");
    }
    if (numParams == 1) {
      cat(hyperparameters[1]);
    } else {
      if (numParams > 4)
        cat("c(", toString(format(hyperparameters[1:4], digits=digits, scientific=TRUE)), ", ...)", sep="")
      else
        cat("c(", toString(format(hyperparameters, digits=digits, scientific=TRUE)), ")", sep="");
    }
    cat(", ", COMMON_SCALE_OPTION_NAME, " = ", commonScaleEnum[commonScale], sep="");
    cat(")");
  } else if (families[1] == getEnumOrder(familyEnum, POINT_FAMILY_NAME)) {
    hyperparameters <- round(hyperparameters, digits);
    location <- hyperparameters;

    cat("(", VALUE_HYPERPARAMETER_NAME, " = ", hyperparameters,
        ", ", POSTERIOR_SCALE_OPTION_NAME, " = ", posteriorScaleEnum[posteriorScale],
        ")", sep = "");
  }
  sink();
  close(stringConnection);


  # families used and the like are only relevant to decompositions with
  # multiple families and a mess of hyperparameters in a single object
  return(list(string = stringResult, numFamiles = numFamiliesUsed,
              numScales = numScalesUsed, numHyperparameters = numHyperparametersUsed));
}
