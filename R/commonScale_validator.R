adjustCommonScalePriorScales <- function(regression, prior)
{
  return(prior);
}

validateCommonScalePrior <- function(regression, prior)
{
  errorPrefix <- "Error applying prior to common scale: ";
  
  if (prior$family == POINT_FAMILY_NAME) {
    if (!is.numeric(prior$value)) stop(errorPrefix, "point prior value must be numeric.");
    if (prior$value <= 0) stop(errorPrefix, "point prior value must be positive.");

    if (prior$posteriorScale != SD_SCALE_NAME &&
        prior$posteriorScale != VAR_SCALE_NAME) {
      stop(errorPrefix, "unrecognized posterior scale '", prior$posteriorScale, "'.");
    }
    
    return(prior);
  } else if (prior$family == INVGAMMA_FAMILY_NAME) {
    if (!is.numeric(prior$shape)) stop(errorPrefix, "inv-gamma prior shape must be numeric.");
    if (prior$shape < 0) stop(errorPrefix, "inv-gamma prior shape must be non-negative.");
    
    if (!is.numeric(prior$scale)) stop(errorPrefix, "inv-gamma prior scale must be numeric.");
    if (prior$scale < 0) stop(errorPrefix, "inv-gamma prior scale must be non-negative.");

    if (prior$posteriorScale != SD_SCALE_NAME &&
        prior$posteriorScale != VAR_SCALE_NAME) {
      stop(errorPrefix, "unrecognized posterior scale '", prior$posteriorScale, "'.");
    }
    
    if (prior$posteriorScale == SD_SCALE_NAME && prior$scale != 0)
      stop(errorPrefix, "inv-gamma prior on sd-scale only conjugate if scale is zero.");

    return(prior);
  }
  
  stop("Internal error: common scale prior of family '",
       prior$family, "' cannot be validated.");
}
