adjustCovariancePriorScales <- function(regression, priors, type,
                                        factorNumber)
{
  if (is.null(priors)) return(NULL);
  
  factorDimension <- nrow(regression@ST[[factorNumber]]);
  factorName <- getFactorNameForNumber(regression, factorNumber);
  
  errorPrefix <- paste("Error applying prior to factor '", factorName, "': ", sep="");
  for (i in 1:length(priors)) {
    prior <- priors[[i]];
    if (is.null(prior)) next;

    if (prior$dataScale == ABSOLUTE_SCALE_NAME) next;

    dataSd <- ifelse(isLinearMixedModel(regression), sd(regression@y), 1);
    inputSd <- getInputSdForFactor(regression, factorNumber);
    sdRatio <- dataSd / inputSd;
    
    if (prior$family == WISHART_FAMILY_NAME) {
      if (type == CORRELATION_TYPE_NAME) {
        prior$scale <- prior$scale * sdRatio;
      } else {
        prior$scale <- prior$scale * sdRatio^2;
      }
    } else if (prior$family == INVWISHART_FAMILY_NAME) {
      if (type == CORRELATION_TYPE_NAME) {
        prior$inverseScale <- prior$inverseScale * sdRatio;
      } else {
        prior$inverseScale <- prior$inverseScale * sdRatio^2;
      }
    } else if (prior$family == GAMMA_FAMILY_NAME) {
      if (prior$posteriorScale == SD_SCALE_NAME) {
        if (type == CORRELATION_TYPE_NAME) {
          prior$rate <- prior$rate / sqrt(sdRatio[i]);
        } else {
          prior$rate <- prior$rate / sdRatio[i];
        }
      } else {
        if (type == CORRELATION_TYPE_NAME) {
          prior$rate <- prior$rate / sdRatio[i];
        } else {
          prior$rate <- prior$rate / sdRatio[i]^2;
        }
      }
    } else if (prior$family == INVGAMMA_FAMILY_NAME) {
      if (prior$posteriorScale == SD_SCALE_NAME) {
        if (type == CORRELATION_TYPE_NAME) {
          prior$scale <- prior$scale * sqrt(sdRatio[i]);
        } else {
          prior$scale <- prior$scale * sdRatio[i];
        }
      } else {
        if (type == CORRELATION_TYPE_NAME) {
          prior$scale <- prior$scale * sdRatio[i];
        } else {
          prior$scale <- prior$scale * sdRatio[i]^2;
        }
      }
    }
    priors[[i]] <- prior;
  }
  return(priors);
}

validateCovariancePrior <- function(regression, priors, factorNumber)
{
  if (is.null(priors)) return(NULL); # occurs for flat priors
  
  factorDimension <- nrow(regression@ST[[factorNumber]]);
  factorName <- getFactorNameForNumber(regression, factorNumber);
  
  errorPrefix <- paste("Error applying prior to factor '", factorName, "': ", sep="");
  for (i in 1:length(priors)) {
    prior <- priors[[i]];
    if (is.null(prior)) next;
    
    if (prior$dataScale != ABSOLUTE_SCALE_NAME &&
        prior$dataScale != FREE_SCALE_NAME) {
      stop(errorPrefix,  "data scale must be '",
           ABSOLUTE_SCALE_NAME, "' or '", FREE_SCALE_NAME, "', was '",
           prior$dataScale, "'.");
    }
    
    if (prior$family == WISHART_FAMILY_NAME) {
      # check degrees of freedom
      degreesOfFreedom <- prior$degreesOfFreedom;
      
      if (!is.numeric(degreesOfFreedom))
        stop(errorPrefix, "degrees of freedom for Wishart evaluated as non-numeric.");
      if (length(degreesOfFreedom) != 1)
        stop(errorPrefix, "degrees of freedom for Wishart evaluated to a numeric value with length != 1.");
      
      if (degreesOfFreedom <= factorDimension - 1)
        stop(errorPrefix, "degrees of freedom for Wishart must be greater the factor dimension - 1.");
      
      
      # check scale matrix
      scale <- prior$scale;
      if (!is.numeric(scale))
        stop(errorPrefix, "scale for Wishart evaluated as non-numeric.");
      
      if (!is.matrix(scale)) {
        if (length(scale) != 1 && length(scale) != factorDimension)
          stop(errorPrefix, "scale for Wishart evaluated to a numeric value with unsuitable length.");
        scale <- diag(scale, factorDimension);
      }
      
      if (nrow(scale) != ncol(scale))
        stop(errorPrefix, "scale for Wishart not a square matrix.");
      
      if (nrow(scale) != factorDimension)
        stop(errorPrefix, "scale for Wishart has dimensions not equal to that of the factor.");
      
      if (!isSymmetric(scale))
        stop(errorPrefix, "scale for Wishart is not symmetric.");
      
      logDetScale <- determinant(scale, TRUE);
      
      if (logDetScale$sign < 0 || logDetScale$modulus == -Inf)
        stop(errorPrefix, "scale for Wishart is not positive definite.");

      priors[[i]] <- list(family = WISHART_FAMILY_NAME,
                          degreesOfFreedom = degreesOfFreedom,
                          scale = scale,
                          dataScale = prior$dataScale);
      
    } else if (prior$family == INVWISHART_FAMILY_NAME) {
      degreesOfFreedom <- prior$degreesOfFreedom;
      
      if (!is.numeric(degreesOfFreedom))
        stop(errorPrefix, "degrees of freedom for Inverse Wishart evaluated as non-numeric.");
       
      if (length(degreesOfFreedom) != 1)
        stop(errorPrefix, "degrees of freedom for Inverse Wishart evaluated to a numeric value with length != 1.");
      
      if (degreesOfFreedom <= factorDimension - 1)
        stop(errorPrefix, "degrees of freedom for Inverse Wishart must be greater than the factor dimension - 1.");
      
      
      inverseScale <- prior$inverseScale;
      
      if (!is.numeric(inverseScale))
        stop(errorPrefix, "inverse scale for Inverse Wishart evaluated as non-numeric.");
      
      if (!is.matrix(inverseScale)) {
        if (length(inverseScale) != 1 && length(inverseScale) != factorDimension)
          stop(errorPrefix, "inverse scale for Inverse Wishart evaluated to a numeric value with unsuitable length.");
        inverseScale <- diag(inverseScale, factorDimension);
      }
      
      if (nrow(inverseScale) != ncol(inverseScale))
        stop(errorPrefix, "inverse scale for Inverse Wishart not a square matrix.");
      
      if (nrow(inverseScale) != factorDimension)
        stop(errorPrefix, "inverse scale for Inverse Wishart has dimensions not equal to that of the factor.");
      
      if (!isSymmetric(inverseScale))
        stop(errorPrefix, "inverse scale for Inverse Wishart is not symmetric.");
      
      logDetInverseScale <- determinant(inverseScale, TRUE);
      
      if (logDetInverseScale$sign < 0 || logDetInverseScale$modulus == -Inf)
        stop(errorPrefix, "inverse scale for Inverse Wishart is not positive definite.");
      
      priors[[i]] <- list(family = INVWISHART_FAMILY_NAME,
                         degreesOfFreedom = degreesOfFreedom,
                         inverseScale = inverseScale,
                          dataScale = prior$dataScale);
    } else if (prior$family == GAMMA_FAMILY_NAME) {
      if (!is.numeric(prior$shape)) stop("gamma shape must be numeric.");
      if (prior$shape <= 0) stop(errorPrefix, "gamma shape must be greater than 0.");
      if (!is.numeric(prior$rate)) stop("gamma rate must be numeric.");
      if (prior$rate <= 0) stop(errorPrefix, "gamma rate must be greater than 0.");
    } else if (prior$family == INVGAMMA_FAMILY_NAME) {
      if (!is.numeric(prior$shape)) stop("inverse Gamma shape must be numeric.");
      if (prior$shape <= 0) stop(errorPrefix, "inverse Gamma shape must be greater than 0.");
      if (!is.numeric(prior$scale)) stop("inverse Gamma scale must be numeric.");
      if (prior$scale <= 0) stop(errorPrefix, "inverse Gamma scale must be greater than 0.");
    }
  }
  
  return(priors);
}
