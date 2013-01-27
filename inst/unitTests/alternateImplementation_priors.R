blme_getPriorPenalty <- function(model) {
  dimpgamma <- function(x, shape, rate, log = TRUE) {
    useLog <- log; rm(log);
    if (shape < 0 || rate < 0) stop("Gamma requires non-negative parameters.");
    if (shape == 0 && rate == 0) stop("Srsly?");
    if (shape == 0) {
      result <- -rate * x;
    } else if (rate == 0) {
      result <- (shape - 1.0) * log(x);
    } else {
      result <- shape * log(rate) + (shape - 1.0) * log(x) - rate * x - lgamma(shape);
    }
    
    return(ifelse(useLog, result, exp(result)));
  }
  dimpigamma <- function(x, shape, scale, log = TRUE) {
    useLog <- log; rm(log);
    if (shape < 0 || scale < 0) stop("Inv-Gamma requires non-negative parameters.");
    if (shape == 0 && scale == 0) stop("Inv-Gamma requires at least one positive parameter.");
    if (shape == 0) {
      result <- -scale / x;
    } else if (scale == 0) {
      result <- -(shape + 1.0) * log(x);
    } else {
      result <- shape * log(scale) - (shape + 1.0) * log(x) - scale / x - lgamma(shape);
    }

    return(ifelse(useLog, result, exp(result)));
  }

  dwish <- function(W, df, logDetScale, scaleInverse, log = FALSE) {
    useLog <- log; rm(log);

    d <- nrow(W);
    
    if (df <  d - 1 && is.infinite(logDetScale)) stop("Wishart requires non-negative parameters.");
    if (df == d - 1 && is.infinite(logDetScale)) stop("Wishart requires at least one positive parameter.");

    result <- -0.5 * sum(W * scaleInverse);
    if (df > d - 1) {
      logDetW <- as.numeric(determinant(W, logarithm=TRUE)$modulus);
      result <- result + 0.5 * (df - d - 1.0) * logDetW;

      if (is.finite(logDetScale)) {
        denominator <- 0.5 * (df * (d * log(2.0) + logDetScale) + 0.5 * d * (d - 1) * log(pi));
        for (i in 1:d) {
          denominator <- denominator + lgamma(0.5 * (df + 1.0 - i));
        }
        result <- result - denominator;
      }
    }
    return(ifelse(useLog, result, exp(result)));
  }
  
  diwish <- function(W, df, logDetInvScale, inverseScale, log=FALSE) {
    useLog <- log; rm(log);

    d <- nrow(W);
        
    if (df < d - 1 && is.infinite(logDetInvScale)) stop("Inv-Wishart requires non-negative parameters.");
    if (df == d - 1 && is.infinite(logDetInvScale)) stop("Inv-Wishart requires at least one positive parameter.");

    result <- 0;
    if (is.finite(logDetInvScale)) {
      W.inv <- solve(W);
    
      result <- result - 0.5 * sum(W.inv * inverseScale);
    }
    if (df > d - 1) {
      logDetW <- as.numeric(determinant(W, logarithm=TRUE)$modulus);
      result <- result - 0.5 * (df + d + 1.0) * logDetW;
    }
    if (df > d - 1 && is.finite(logDetInvScale)) {
      denominator <- 0.5 * (df * (d * log(2.0) - logDetInvScale) + 0.5 * (d - 1) * log(pi));
      for (i in 1:d) {
        denominator <- denominator + lgamma(0.5 * (df + 1.0 - i));
      }
      result <- result - denominator;
    }
    return(ifelse(useLog, result, exp(result)));
  }
    
  
  result <- 0;
  commonScale <- ifelse(model@dims[["REML"]] == 1, model@deviance[["sigmaREML"]], model@deviance[["sigmaML"]])^2;
    

  directType <- blme:::getEnumOrder(blme:::typeEnum, blme:::DIRECT_TYPE_NAME);
  
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


  if (model@var.prior@type  == directType) {
    prior <- model@var.prior;
    family <- prior@families[1];
    if (family == gammaFamily) {
      shape <- prior@hyperparameters[1];
      rate  <- prior@hyperparameters[2];
      onVarScale <- any(prior@scales[1] == c(varOnCommonScale, varNotOnCommonScale));
      if (onVarScale) {
        result <- result - 2.0 * dimpgamma(commonScale, shape, rate, TRUE);
      } else {
        result <- result - 2.0 * dimpgamma(sqrt(commonScale), shape, rate, TRUE);
      }
    } else if (family == invGammaFamily) {
      shape <- prior@hyperparameters[1];
      scale <- prior@hyperparameters[2];
      onVarScale <- any(prior@scales[1] == c(varOnCommonScale, varNotOnCommonScale));
      if (onVarScale) {
        result <- result - 2.0 * dimpigamma(commonScale, shape, scale, TRUE);
      } else {
        result <- result - 2.0 * dimpigamma(sqrt(commonScale), shape, scale, TRUE);
      }
    }
  }

  if (model@fixef.prior@type == directType) {
    prior <- model@fixef.prior;
    family <- prior@families[1];
    numUnmodeledCoef <- model@dims[["p"]];
    
    if (family == normalFamily) {
      exponentialTerm <- 0;
      useCommonScale <- any(prior@scales[1] == c(varOnCommonScale, sdOnCommonScale));

      logDetCov <- prior@hyperparameters[1];
      numHyperparameters <- length(prior@hyperparameters) - 1;
      hyperparameters <- prior@hyperparameters[2:length(prior@hyperparameters)];
      
      if (numHyperparameters == 1) {
        varInverse <- hyperparameters[1]^2;
        exponentialTerm <- crossprod(model@fixef)[1] * varInverse;
      } else if (numHyperparameters == numUnmodeledCoef) {
        varsInverse <- hyperparameters^2;
        
        for (i in 1:numUnmodeledCoef) {
          exponentialTerm <- exponentialTerm + model@fixef[i]^2 * varsInverse[i];
        }
      } else {
        covLength <- numUnmodeledCoef^2;
        covInverse <- matrix(hyperparameters[covLength + 1:covLength], numUnmodeledCoef, numUnmodeledCoef);
        exponentialTerm <- (crossprod(model@fixef, covInverse) %*% model@fixef)[1];
      }

      if (useCommonScale) {
        exponentialTerm <- exponentialTerm / commonScale;
        logDetCov <- logDetCov + numUnmodeledCoef * log(commonScale);
      }
      
      result <- result + numUnmodeledCoef * log(2 * pi) + logDetCov + exponentialTerm;
    }
  }
    
  for (i in 1:length(model@ST)) {
    prior <- model@cov.prior[[i]];
    if (prior@type != directType) next; # have not implemented correlation or spectral decomps
    
    family <- prior@families[1];
    
    covar <- blme_stToCov(model@ST[[i]]);
    if (family == gammaFamily) {
      shape <- prior@hyperparameters[1];
      rate  <- prior@hyperparameters[2];
      param <- covar[1];
      
      onCommonScale <- any(prior@scales[1] == c(varOnCommonScale, sdOnCommonScale));
      onVarScale    <- any(prior@scales[1] == c(varOnCommonScale, varNotOnCommonScale));
      
      if (!onCommonScale) param <- param * commonScale;
      if (!onVarScale) param <- sqrt(param);
      
      result <- result - 2.0 * dimpgamma(param, shape, rate, TRUE);
    } else if (family == invGammaFamily) {
      shape <- prior@hyperparameters[1];
      scale <- prior@hyperparameters[2];
      param <- covar[1];
      
      onCommonScale <- any(prior@scales[1] == c(varOnCommonScale, sdOnCommonScale));
      onVarScale    <- any(prior@scales[1] == c(varOnCommonScale, varNotOnCommonScale));
      
      if (!onCommonScale) param <- param * commonScale;
      if (!onVarScale) param <- sqrt(param);
      
      result <- result - 2.0 * dimpigamma(param, shape, scale, TRUE);
    } else if (family == wishartFamily) {
      factorDimension <- nrow(covar);
      df <- prior@hyperparameters[1];
      logDetScale <- prior@hyperparameters[2];
      scaleInverse <- matrix(prior@hyperparameters[-c(1,2)], factorDimension, factorDimension);
      
      onCommonScale <- any(prior@scales[1] == c(varOnCommonScale, sdOnCommonScale));
      if (!onCommonScale) covar <- covar * commonScale;
      
      result <- result - 2.0 * dwish(covar, df, logDetScale, scaleInverse, TRUE);
    } else if (family == invWishartFamily) {
      factorDimension <- nrow(covar);
      df <- prior@hyperparameters[1];
      logDetInverseScale <- prior@hyperparameters[2];
      inverseScale <- matrix(prior@hyperparameters[-c(1,2)], factorDimension, factorDimension)
      
      onCommonScale <- any(prior@scales[1] == c(varOnCommonScale, sdOnCommonScale));
      if (!onCommonScale) covar <- covar * commonScale;
      
      result <- result - 2.0 * diwish(covar, df, logDetInverseScale, inverseScale, TRUE);
    }
  }
  return(result);
}
