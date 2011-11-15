adjustCommonScalePriorScales <- function(regression, prior)
{
  return(prior);
}

validateCommonScalePrior <- function(regression, prior)
{
  errorPrefix <- paste("Error applying prior to common scale: ", sep="");
  
  if (prior$family == POINT_FAMILY_NAME) {
    if (!is.numeric(prior$value)) stop(errorPrefix, "point prior value must be numeric.");
    if (prior$value <= 0) stop(errorPrefix, "point prior value must be positive.");
    
    return(prior);
  }
  
  stop("Internal error: common scale prior of family '",
       prior$family, "' cannot be validated.");
}
