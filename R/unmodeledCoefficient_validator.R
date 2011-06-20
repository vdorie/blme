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
      prior$dataScale != FREE_SCALE_NAME) {
    stop(errorPrefix,  "data scale must be '",
         ABSOLUTE_SCALE_NAME, "' or '", FREE_SCALE_NAME, "', was '",
         prior$dataScale, "'.");
  }
  
  if (prior$family == NORMAL_FAMILY_NAME) {
    covarianceScale <- prior[["covarianceScale"]];
    if (covarianceScale != ABSOLUTE_SCALE_NAME &&
        covarianceScale != COMMON_SCALE_NAME)
      stop(errorPrefix, "covariance scale must be '",
         ABSOLUTE_SCALE_NAME, "' or '", COMMON_SCALE_NAME, "', was '",
         covarianceScale, "'.");

    if (!isLinearMixedModel(regression) && covarianceScale == COMMON_SCALE_NAME)
      stop(errorPrefix, "covariance scale of '", COMMON_SCALE_NAME,
           "' can only be applied to linear mixed models.");

    
    covariance <- prior[["covariance"]];

    if (!is.numeric(covariance)) stop(errorPrefix, "covariance for Normal evaluated as non-numeric.");

    if (!is.matrix(covariance)) {
      if (length(covariance) != 1 && length(covariance) != 2 &&
          length(scale) != numUnmodeledCoefs)
        stop(errorPrefix, "scale for Wishart evaluated to a numeric value with unsuitable length.");

      if (any(covariance <= 0)) stop(errorPrefix, "covariance for Normal is not positive definite.");
      
      if (length(covariance) == 2) {
        covariance <- diag(c(covariance[1], rep(covariance[2], numUnmodeledCoefs - 1)), numUnmodeledCoefs);
        prior$covariance <- covariance;
      }
    }
    
    if (is.matrix(covariance)) {  
      if (nrow(covariance) != ncol(covariance))
        stop(errorPrefix, "covariance for Normal not a square matrix.");
      
      if (nrow(covariance) != numUnmodeledCoefs)
        stop(errorPrefix, "covariance for Normal has dimensions not equal to the number of unmodeled coefficients.");
      
      if (!isSymmetric(covariance))
        stop(errorPrefix, "covariance for Normal is not symmetric.");
      
      logDetCov <- determinant(covariance, TRUE);
      
      if (logDetCov$sign < 0 || logDetCov$modulus == -Inf)
        stop(errorPrefix, "covariance for Normal is not positive definite.");
    }   
    
    return(prior);
  }
  
  stop("Internal error: unmodeled coefficient prior of family '",
       prior$family, "' cannot be validated.");
}
