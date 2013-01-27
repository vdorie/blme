blme_getPriorDegreesOfFreedom <- function(model) {
  noneType    <- blme:::getEnumOrder(blme:::typeEnum, blme:::NONE_TYPE_NAME);
  directType  <- blme:::getEnumOrder(blme:::typeEnum, blme:::DIRECT_TYPE_NAME);
  #corrType <- blme:::getEnumOrder(blme:::familyEnum, blme:::CORRELATION_TYPE_NAME);
  #specType <- blme:::getEnumOrder(blme:::familyEnum, blme:::SPECTRAL_TYPE_NAME);
  
  gammaFamily <- blme:::getEnumOrder(blme:::familyEnum, blme:::GAMMA_FAMILY_NAME);
  invGammaFamily <- blme:::getEnumOrder(blme:::familyEnum, blme:::INVGAMMA_FAMILY_NAME);
  wishartFamily    <- blme:::getEnumOrder(blme:::familyEnum, blme:::WISHART_FAMILY_NAME);
  invWishartFamily <- blme:::getEnumOrder(blme:::familyEnum, blme:::INVWISHART_FAMILY_NAME);
  normalFamily <- blme:::getEnumOrder(blme:::familyEnum, blme:::NORMAL_FAMILY_NAME);
  
  varOnCommonScale <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::VAR_SCALE_NAME),
                                         blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_TRUE_NAME));
  sdOnCommonScale  <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::SD_SCALE_NAME),
                                         blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_TRUE_NAME));
  varNotOnCommonScale <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::VAR_SCALE_NAME),
                                            blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_FALSE_NAME));
  
  getCommonScalePriorDF <- function(model) {
    if (model@var.prior@type[1] != directType) return(0.0);
    
    family <- model@var.prior@families[1];
    hyperparameters <- model@var.prior@hyperparameters;
    
    # obviously should be true, but to be consistent
    onCommonScale <- blme:::getEnumOrder(blme:::commonScaleEnum, blme:::defaultCommonScaleCommonScale);
    onVarScale    <- blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::VAR_SCALE_NAME);
    
    onVarScale <- model@var.prior@scales[1] == blme:::getScaleInt(onVarScale, onCommonScale);

    result <- 0;
    if (family == gammaFamily) {
      shape <- hyperparameters[1];
        
      result <- result - ifelse(onVarScale, 2.0 * (shape - 1.0), shape - 1.0);
    } else if (family == invGammaFamily) {
      shape <- hyperparameters[1];
      
      result <- result + ifelse(onVarScale, 2.0 * (shape + 1.0), shape + 1.0);
    }
    return(result);
  }

  getFixefPriorDF <- function(model) {
    result <- 0;
    if (model@fixef.prior@type[1]     == directType &&
        model@fixef.prior@families[1] == normalFamily) {
      onVarScale    <- blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::defaultUnmodeledCoefficientPosteriorScale);
      onCommonScale <- blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_TRUE_NAME);
      
      onCommonScale <- model@fixef.prior@scales[1] == blme:::getScaleInt(onVarScale, onCommonScale);
      
      if (onCommonScale) result <- result + model@dims[["p"]];
    }
    return(result);
  }

  getCovPriorDF <- function(model) {
    getDFForCovPrior <- function(family, scale, hyperparameters, factorDimension) {
      onCommonScale <- scale == varOnCommonScale || scale == sdOnCommonScale;
      onVarScale    <- scale == varOnCommonScale || scale == varNotOnCommonScale;
      
      result <- 0.0;
      numHyperparamsUsed <- 0;
      if (family == gammaFamily) {
        shape <- hyperparameters[1];
        numHyperparametersUsed <- 2;
        
        result <- ifelse(onVarScale, -2.0 * (shape - 1.0), -(shape - 1.0));
      } else if (family == invGammaFamily) {
        shape <- hyperparameters[1];
        numHyperparametersUsed <- 2;
        
        result <- ifelse(onVarScale, 2.0 * (shape + 1.0), shape + 1.0);
      } else if (family == wishartFamily) {
        df <- hyperparameters[1];
        numHyperparametersUsed <- factorDimension * factorDimension;
        
        result <- -factorDimension * (df - factorDimension - 1.0);
      } else if (family == invWishartFamily) {
        df <- hyperparameters[1];
        numHyperparametersUsed <- factorDimension * factorDimension;
        
        result <- factorDimension * (df + factorDimension + 1.0);
      }
      if (onCommonScale) return(list(df = 0.0, numHyperparametersUsed = numHyperparametersUsed));
      return(list(df = result, numHyperparametersUsed = numHyperparametersUsed));
    }
    
    result <- 0.0;
    for (i in 1:model@dims[["nt"]]) {
      factorDimension <- nrow(model@ST[[i]]);
      prior <- model@cov.prior[[i]];

      if (prior@type == noneType) next;
      
      loopMax <- ifelse(prior@type == directType, 1, factorDimension);
      hyperparameters <- prior@hyperparameters;
      for (j in 1:loopMax) {
        family <- prior@families[j];
        scale  <- prior@scales[j];

        priorDfs <- getDFForCovPrior(family, scale, hyperparameters, factorDimension);
        result <- result + priorDfs$df;
        hyperparameters <- hyperparameters[-(1:priorDfs$numHyperparametersUsed)];
      }
    }
    return(result);
  }
  
  return(getCommonScalePriorDF(model) + getFixefPriorDF(model) + getCovPriorDF(model));  
}


blme_getCommonScaleExponentialPart <- function(model, parameters) {
  noneType   <- blme:::getEnumOrder(blme:::typeEnum, blme:::NONE_TYPE_NAME);
  directType <- blme:::getEnumOrder(blme:::typeEnum, blme:::DIRECT_TYPE_NAME);
  corrType   <- blme:::getEnumOrder(blme:::typeEnum, blme:::CORRELATION_TYPE_NAME);
  specType   <- blme:::getEnumOrder(blme:::typeEnum, blme:::SPECTRAL_TYPE_NAME);
  
  gammaFamily      <- blme:::getEnumOrder(blme:::familyEnum, blme:::GAMMA_FAMILY_NAME);
  invGammaFamily   <- blme:::getEnumOrder(blme:::familyEnum, blme:::INVGAMMA_FAMILY_NAME);
  wishartFamily    <- blme:::getEnumOrder(blme:::familyEnum, blme:::WISHART_FAMILY_NAME);
  invWishartFamily <- blme:::getEnumOrder(blme:::familyEnum, blme:::INVWISHART_FAMILY_NAME);
  normalFamily     <- blme:::getEnumOrder(blme:::familyEnum, blme:::NORMAL_FAMILY_NAME);
  
  varOnCommonScale <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::VAR_SCALE_NAME),
                                         blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_TRUE_NAME));
  sdOnCommonScale  <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::SD_SCALE_NAME),
                                         blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_TRUE_NAME));
  varNotOnCommonScale <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::VAR_SCALE_NAME),
                                            blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_FALSE_NAME));
  
  getCommonScaleExponentialPart <- function(model) {
    result <- list(mTwo = 0, mOne = 0, two = 0, one = 0);

    if (model@var.prior@type[1] != directType) return(result);
     
    family <- model@var.prior@families[1];
    hyperparameters <- model@var.prior@hyperparameters;
    
    # obviously should be true, but to be consistent
    onCommonScale <- blme:::getEnumOrder(blme:::commonScaleEnum, blme:::defaultCommonScaleCommonScale);
    onVarScale    <- blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::VAR_SCALE_NAME);
    
    onVarScale <- model@var.prior@scales[1] == blme:::getScaleInt(onVarScale, onCommonScale);
      
    if (family == gammaFamily) {
      rate <- hyperparameters[2];
      if (onVarScale) {
        result$two <- rate;
      } else {
        result$one <- rate;
      }
    } else if (family == invGammaFamily) {
      scale <- hyperparameters[2];
      if (onVarScale) {
        result$mTwo <- scale;
      } else {
        result$mOne <- scale;
      }
    }
    return(result);
  }

  getCovarianceExponentialPart <- function(model, parameters) {
    # assumes parameters are on sd scale for univariate
    getExponentialPartForCovPrior <- function(family, scale, hyperparameters, parameters, factorDimension) {
      result <- list(mTwo = 0, mOne = 0, two = 0, one = 0,
                     numHyperparametersUsed = 0, numParametersUsed = 0);
      
      onCommonScale <- scale == varOnCommonScale || scale == sdOnCommonScale;
      onVarScale    <- scale == varOnCommonScale || scale == varNotOnCommonScale;
      numHyperparamsUsed <- 0;
      if (family == gammaFamily) {
        rate <- hyperparameters[2];
        result$numHyperparametersUsed <- 2; result$numParametersUsed <-  1;

        if (onCommonScale) return(result);
        if (onVarScale) {
          result$two <- result$two + rate * parameters[1]^2;
        } else {
          result$one <- result$one + rate * parameters[1];
        }
      } else if (family == invGammaFamily) {
        scale <- hyperparameters[1];
        result$numHyperparametersUsed <- 2; result$numParametersUsed <-  1;
        
        if (onCommonScale) return(result);
        if (onVarScale) {
          result$mTwo <- result$mTwo + scale / parameters[1]^2;
        } else {
          result$mOne <- result$mOne + scale / parameters[1];
        }
      } else if (family == wishartFamily) {
        scaleInverse <- matrix(hyperparameters[1:factorDimension^2 + 2], factorDimension, factorDimension);
        result$numHyperparametersUsed <- 2 + factorDimension^2; result$numParametersUsed <- result$numParametersUsed + factorDimension * (factorDimension + 1) / 2;

        if (onCommonScale) return(result);
        
        cov <- blme_stToCov(blme_stVectorToMatrix(parameters, factorDimension));
        result$two <- result$two + 0.5 * sum(cov * scaleInverse); # trace of crossproduct of symmetric matrices
      } else if (family == invWishartFamily) {
        inverseScale <- matrix(hyperparameters[1:factorDimension^2 + 2], factorDimension, factorDimension);
        result$numHyperparametersUsed <- 2 + factorDimension^2; result$numParametersUsed <- result$numParametersUsed + factorDimension * (factorDimension + 1) / 2;

        if (onCommonScale) return(result);

        covInverse <- blme_stToCovInv(blme_stVectorToMatrix(parameters, factorDimension));
        result$mTwo <- result$mTwo + 0.5 * sum(covInverse * inverseScale);
      }
      return(result);
    }
    
    result <- list(mTwo = 0, mOne = 0, two = 0, one = 0);
    for (i in 1:model@dims[["nt"]]) {
      factorDimension <- nrow(model@ST[[i]]);
      prior <- model@cov.prior[[i]];

      if (prior@type == noneType) next;
      
      loopMax <- ifelse(prior@type == directType, 1,
                        ifelse(prior@type == specType, factorDimension,
                               factorDimension + 1));
      hyperparameters <- prior@hyperparameters;
      for (j in 1:loopMax) {
        family <- prior@families[j];
        scale  <- prior@scales[j];

        priorExpPart <- getExponentialPartForCovPrior(family, scale, hyperparameters, parameters, factorDimension);
        for (name in names(result)) result[[name]] <- result[[name]] + priorExpPart[[name]];
        
        hyperparameters <- hyperparameters[-(1:priorExpPart$numHyperparametersUsed)];
        parameters <- parameters[-(1:priorExpPart$numParametersUsed)];
      }
    }
    return (result);
  }
  
  result <- getCommonScaleExponentialPart(model);
  covResult <- getCovarianceExponentialPart(model, parameters);
  for (name in names(result)) result[[name]] <- result[[name]] + covResult[[name]];
  return(result);
}


blme_canProfileCommonScale <- function(model) {
  noneType   <- blme:::getEnumOrder(blme:::typeEnum, blme:::NONE_TYPE_NAME);
  directType <- blme:::getEnumOrder(blme:::typeEnum, blme:::DIRECT_TYPE_NAME);
  corrType   <- blme:::getEnumOrder(blme:::typeEnum, blme:::CORRELATION_TYPE_NAME);
  specType   <- blme:::getEnumOrder(blme:::typeEnum, blme:::SPECTRAL_TYPE_NAME);
  
  gammaFamily      <- blme:::getEnumOrder(blme:::familyEnum, blme:::GAMMA_FAMILY_NAME);
  invGammaFamily   <- blme:::getEnumOrder(blme:::familyEnum, blme:::INVGAMMA_FAMILY_NAME);
  wishartFamily    <- blme:::getEnumOrder(blme:::familyEnum, blme:::WISHART_FAMILY_NAME);
  invWishartFamily <- blme:::getEnumOrder(blme:::familyEnum, blme:::INVWISHART_FAMILY_NAME);
  normalFamily     <- blme:::getEnumOrder(blme:::familyEnum, blme:::NORMAL_FAMILY_NAME);
  pointFamily      <- blme:::getEnumOrder(blme:::familyEnum, blme:::POINT_FAMILY_NAME);
  
  varOnCommonScale <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::VAR_SCALE_NAME),
                                         blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_TRUE_NAME));
  sdOnCommonScale  <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::SD_SCALE_NAME),
                                         blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_TRUE_NAME));
  varNotOnCommonScale <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::VAR_SCALE_NAME),
                                            blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_FALSE_NAME));

  covPriorComplicatesCommonScale <- function(family, scale, hyperparameters, factorDimension) {
    onCommonScale <- scale == varOnCommonScale || scale == sdOnCommonScale;
    onVarScale    <- scale == varOnCommonScale || scale == varNotOnCommonScale;

    result <- FALSE;
    numHyperparamsUsed <- 0;
    if (family == gammaFamily) {
      rate <- hyperparameters[2];
      numHyperparametersUsed <- 2;

      result <- rate > 0.0;
    } else if (family == invGammaFamily) {
      scale <- hyperparameters[2];
      numHyperparametersUsed <- 2;
      
      result <- !onVarScale && scale > 0;
    } else if (family == wishartFamily) {
      logDetScale <- hyperparameters[2];
      numHyperparametersUsed <- factorDimension * factorDimension;

      result <- is.finite(logDetScale);
    } else if (family == invWishartFamily) {
      numHyperparametersUsed <- factorDimension * factorDimension;

      result <- FALSE;
    }
    if (onCommonScale) return(list(result = FALSE, numHyperparametersUsed = numHyperparametersUsed));
    return(list(result = result, numHyperparametersUsed = numHyperparametersUsed));
  }
  
  commonScalePrior <- model@var.prior;
  if (commonScalePrior@type == directType) {
    family <- commonScalePrior@families[1];
    if (family == pointFamily) { return(FALSE); }
    else if (family == invGammaFamily) {
      scale <- commonScalePrior@hyperparameters[2];
      postScale <- commonScalePrior@scales[1];
      if (postScale != varOnCommonScale && postScale != varNotOnCommonScale) {
        # sd prior, on profile if scale is 0
        if (scale > 0) return(FALSE);
      }
    } else if (family == gammaFamily) {
      rate <-  commonScalePrior@hyperparameters[2];
      if (rate > 0) return(FALSE);
    }
  }

  fixefPrior <- model@fixef.prior;
  if (fixefPrior@type == directType) {
    fixefPriorIsNormal <- fixefPrior@families[1] == normalFamily;
    
    onCommonScale <- FALSE;
    if (fixefPriorIsNormal) {
      onCommonScale <- any(fixefPrior@scales[1] == c(varOnCommonScale, sdOnCommonScale));
    }
    if (fixefPriorIsNormal && !onCommonScale) return(FALSE);
  }

  for (i in 1:model@dims[["nt"]]) {
    factorDimension <- nrow(model@ST[[i]]);
    prior <- model@cov.prior[[i]];
    
    if (prior@type == noneType) next;
    
    loopMax <- ifelse(prior@type == directType, 1, factorDimension);
    hyperparameters <- prior@hyperparameters;
    for (j in 1:loopMax) {
      family <- prior@families[j];
      scale  <- prior@scales[j];

      familyResult <- covPriorComplicatesCommonScale(family, scale, hyperparameters, factorDimension);
      if (familyResult$result == TRUE) return(FALSE);
      hyperparameters <- hyperparameters[-(1:familyResult$numHyperparametersUsed)];
    }
  }
  return(TRUE);
}

blme_parametersIncludeCommonScale <- function(model)
{
  if (blme_canProfileCommonScale(model)) return(FALSE);

  directType  <- blme:::getEnumOrder(blme:::typeEnum, blme:::DIRECT_TYPE_NAME);
  pointFamily <- blme:::getEnumOrder(blme:::familyEnum, blme:::POINT_FAMILY_NAME);
  
  commonScalePrior <- model@var.prior;
  if (commonScalePrior@type == directType) {
    family <- commonScalePrior@families[1];
    if (family == pointFamily) { return(FALSE); }
  }

  return(TRUE);
}
