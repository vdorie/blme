# functions that the prior parser might call that don't influence the logic of the
# of the parse
getEnumOrder <- function(enumeration, name)
{
  return (as.integer(which(enumeration == name) - 1));
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
      if (result[i] == 0) result[i] <- 1;
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


buildStringForFamily <- function(families, scales, hyperparameters, preprocessed)
{
  matchedCall <- match.call();
  
  numFamiliesUsed <- 0;
  numScalesUsed <- 0;
  numHyperparametersUsed <- 0;

  stringResult <- NULL; # for R CMD Check
  
  stringConnection <- textConnection("stringResult", "w", local=TRUE);
  sink(stringConnection);

  cat(familyEnumeration[families[1] + 1]);

  if (families[1] == getEnumOrder(familyEnumeration, GAMMA_FAMILY_NAME)) {
    cat("(", SHAPE_HYPERPARAMETER_NAME, " = ", hyperparameters[1],
        ", ", RATE_HYPERPARAMETER_NAME, " = ", hyperparameters[2],
        ", ", POSTERIOR_SCALE_OPTION_NAME, " = ", scaleEnumeration[scales[1] + 1],
        ")", sep="");
    numFamiliesUsed <- 1;
    numScalesUsed <- 1;
    numHyperparametersUsed <- 2;
  } else if (families[1] == getEnumOrder(familyEnumeration, INVGAMMA_FAMILY_NAME)) {
    cat("(", SHAPE_HYPERPARAMETER_NAME, " = ", hyperparameters[1],
        ", ", SCALE_HYPERPARAMETER_NAME, " = ", hyperparameters[2],
        ", ", POSTERIOR_SCALE_OPTION_NAME, " = ", scaleEnumeration[scales[1] + 1],
        ")", sep="");
    numFamiliesUsed <- 1;
    numScalesUsed <- 1;
    numHyperparametersUsed <- 2;
  } else if (families[1] == getEnumOrder(familyEnumeration, WISHART_FAMILY_NAME)) {
    scaleInverse <- hyperparameters[3:length(hyperparameters)];
    levelDimension <- round(sqrt(length(scaleInverse)), digits=0);
    scaleInverse <- matrix(scaleInverse, levelDimension, levelDimension);
    scale <- as.numeric(solve(scaleInverse));
    
    
    cat("(", DEGREES_OF_FREEDOM_HYPERPARAMETER_NAME, " = ", hyperparameters[1], ", ",
        SCALE_HYPERPARAMETER_NAME, " = ", sep="");
    if (levelDimension == 1) {
      cat(scale[1]);
    } else {
      if (levelDimension > 2)
        cat("c(", toString(format(scale[1:4], digits=2, scientific=TRUE)), ", ...)", sep="")
      else
        cat("c(", toString(format(scale, digits=2, scientific=TRUE)), ")", sep="");
    }
    cat(")");
  } else if (families[1] == getEnumOrder(familyEnumeration, INVWISHART_FAMILY_NAME)) {
    inverseScale <- hyperparameters[3:length(hyperparameters)];
    levelDimension <- round(sqrt(length(inverseScale)), digits=0);
    
    cat("(", DEGREES_OF_FREEDOM_HYPERPARAMETER_NAME, " = ", hyperparameters[1], ", ",
        INVERSE_SCALE_HYPERPARAMETER_NAME, " = ", sep="");
    if (levelDimension == 1) {
      cat(inverseScale[1]);
    } else {
      if (levelDimension > 2)
        cat("c(", toString(format(inverseScale[1:4], digits=2, scientific=TRUE)), ", ...)", sep="")
      else
        cat("c(", toString(format(inverseScale, digits=2, scientific=TRUE)), ")", sep="");
    }
    cat(")");
  } else if (families[1] == getEnumOrder(familyEnumeration, NORMAL_FAMILY_NAME)) {
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
        cat("c(", toString(format(hyperparameters[1:4], digits=2, scientific=TRUE)), ", ...)", sep="")
      else
        cat("c(", toString(format(hyperparameters, digits=2, scientific=TRUE)), ")", sep="");
    }
    cat(")");
  } else if (families[1] == getEnumOrder(familyEnumeration, POINT_FAMILY_NAME)) {
    location <- hyperparameters;

    cat("(", VALUE_HYPERPARAMETER_NAME, " = ", hyperparameters,
        ", ", POSTERIOR_SCALE_OPTION_NAME, " = ", scaleEnumeration[scales + 1],
        ")", sep = "");
  }
  sink();
  close(stringConnection);


  # families used and the like are only relevant to decompositions with
  # multiple families and a mess of hyperparameters in a single object
  return(list(string = stringResult, numFamiles = numFamiliesUsed,
              numScales = numScalesUsed, numHyperparameters = numHyperparametersUsed));
}
