adjustUnmodeledCoefficientPriorScales <- function(regression, prior)
{
  if (prior$dataScale == ABSOLUTE_SCALE_NAME) return(prior);

  numUnmodeledCoefs <- regression@dims[["p"]];

  dataSd <- ifelse(isLinearMixedModel(regression), sd(regression@y), 1);
  inputSd <- getInputSdForUnmodeledCoefficients(regression);
  sdRatio <- dataSd / inputSd;

  errorPrefix <- paste("Error applying prior to unmodeled coefficients: ", sep="");
  if (prior$family == NORMAL_FAMILY_NAME) {
    covariance <- prior[["covariance"]];
    if (length(covariance) == 1 || length(covariance) == numUnmodeledCoefs) {
      covariance <- covariance * sdRatio^2;
    } else {
      rescalingMatrix <- diag(sdRatio, numUnmodeledCoefs);
      covariance <- rescalingMatrix %*% covariance %*% rescalingMatrix;
    }
    prior[["covariance"]] <- covariance;
  }
  return(prior);
}

validateUnmodeledCoefficientPrior <- function(regression, prior)
{
  numUnmodeledCoefs <- regression@dims[["p"]];
  
  errorPrefix <- paste("Error applying prior to unmodeled coefficients: ", sep="");
    
  if (prior$dataScale != ABSOLUTE_SCALE_NAME &&
      prior$dataScale != FREE_SCALE_NAME)
  {
    stop(errorPrefix,  "data scale must be '",
         ABSOLUTE_SCALE_NAME, "' or '", FREE_SCALE_NAME, "', was '",
         prior$dataScale, "'.");
  }
  
  if (prior$family == NORMAL_FAMILY_NAME) {
    if (!is.null(prior$onCommonScale)) {
      onCommonScale <- prior$onCommonScale;
    } else {
      covScale <- prior[["covarianceScale"]];
      warning("Option 'cov.scale' for fixef priors has been deprecated. Use 'common.scale' = '",
              COMMON_SCALE_TRUE_NAME, "' or '", COMMON_SCALE_FALSE_NAME, "' instead.");
      if (covScale != ABSOLUTE_SCALE_NAME && covScale != COMMON_SCALE_NAME)
        stop(errorPrefix, "cov.scale must be '",
             ABSOLUTE_SCALE_NAME, "' or '", COMMON_SCALE_NAME, "', was '",
             covScale, "'.");
      onCommonScale <- ifelse(covScale == COMMON_SCALE_NAME,
                              COMMON_SCALE_TRUE_NAME,
                              COMMON_SCALE_FALSE_NAME);
    }

    if (is.logical(onCommonScale))
      onCommonScale <- ifelse(onCommonScale, COMMON_SCALE_TRUE_NAME, COMMON_SCALE_FALSE_NAME);
    if (onCommonScale != COMMON_SCALE_TRUE_NAME &&
        onCommonScale != COMMON_SCALE_FALSE_NAME) {
      stop(errorPrefix, "common.scale must be '", COMMON_SCALE_TRUE_NAME, "' or '",
           COMMON_SCALE_FALSE_NAME, "'.");
    }
    prior$onCommonScale <- onCommonScale;
    
    covariance <- prior[["covariance"]];
    if (!is.numeric(covariance)) stop(errorPrefix, "covariance for normal evaluated as non-numeric.");
    
    if (!is.matrix(covariance)) {
      if (length(covariance) != 1 && length(covariance) != 2 &&
          length(scale) != numUnmodeledCoefs)
        stop(errorPrefix, "covariance for normal evaluated to a numeric value with unsuitable length.");
      
      if (any(covariance <= 0)) stop(errorPrefix, "covariance for normal is not positive definite.");
      
      if (length(covariance) == 2) {
        covariance <- diag(c(covariance[1], rep(covariance[2], numUnmodeledCoefs - 1)), numUnmodeledCoefs);
        prior$covariance <- covariance;
      }
    }
    
    if (is.matrix(covariance)) {  
      if (nrow(covariance) != ncol(covariance))
        stop(errorPrefix, "covariance for normal not a square matrix.");
      
      if (nrow(covariance) != numUnmodeledCoefs)
        stop(errorPrefix, "covariance for normal has dimensions not equal to the number of unmodeled coefficients.");
      
      if (!isSymmetric(covariance))
        stop(errorPrefix, "covariance for normal is not symmetric.");
      
      logDetCov <- determinant(covariance, TRUE);
      if (logDetCov$sign < 0 || !is.finite(logDetCov$modulus))
        stop(errorPrefix, "covariance for normal is not positive definite and finite.p");
    }
    
    return(prior);
  }
  
  stop("Internal error: unmodeled coefficient prior of family '",
       prior$family, "' cannot be validated.");
}
