loadConstants <- function(env) {
  env$typeEnum   <- .Call(bmer_getTypeEnumeration);
  env$familyEnum <- .Call(bmer_getFamilyEnumeration);
  env$posteriorScaleEnum  <- .Call(bmer_getPosteriorScaleEnumeration);
  env$commonScaleEnum     <- .Call(bmer_getCommonScaleEnumeration);

  # There's no good way to do constants/enums in R (that I know).
  # This reflects how it's done in C but, I'm open to suggestions.
  #
  # WARNING: should match enums in blmer.c. If the order changes there,
  # it has to change here as well
  
  # "types" relate to decompositions
  env$NONE_TYPE_NAME        <- typeEnum[1];
  env$DIRECT_TYPE_NAME      <- typeEnum[2];
  env$CORRELATION_TYPE_NAME <- typeEnum[3];
  env$SPECTRAL_TYPE_NAME    <- typeEnum[4];

  # families of distributions over parameters
  env$FLAT_FAMILY_NAME       <- familyEnum[1]; # for variances
  env$GAMMA_FAMILY_NAME      <- familyEnum[2];
  env$INVGAMMA_FAMILY_NAME   <- familyEnum[3];
  env$WISHART_FAMILY_NAME    <- familyEnum[4];
  env$INVWISHART_FAMILY_NAME <- familyEnum[5];
  
  env$NORMAL_FAMILY_NAME     <- familyEnum[6]; # for fixed effects
  env$MVT_FAMILY_NAME        <- familyEnum[7];

  env$POINT_FAMILY_NAME      <- familyEnum[8]; # for common scale

  # scales determine how the prior shold be applied
  env$SD_SCALE_NAME       <- posteriorScaleEnum[1]; # variances
  env$VAR_SCALE_NAME      <- posteriorScaleEnum[2];

  env$COMMON_SCALE_FALSE_NAME <- commonScaleEnum[1];
  env$COMMON_SCALE_TRUE_NAME  <- commonScaleEnum[2];

  # similar, but not actually passed down to C
  env$FREE_SCALE_NAME     <- "free";
  
  # deprecated
  env$ABSOLUTE_SCALE_NAME <- "absolute"
  env$COMMON_SCALE_NAME   <- "common"
  
  env$SUPPORTED_COVARIANCE_TYPE_NAMES <- c(NONE_TYPE_NAME, CORRELATION_TYPE_NAME, SPECTRAL_TYPE_NAME,
                                           FLAT_FAMILY_NAME, GAMMA_FAMILY_NAME, INVGAMMA_FAMILY_NAME,
                                           WISHART_FAMILY_NAME, INVWISHART_FAMILY_NAME);
  env$SUPPORTED_UNMODELED_COEFFICIENT_FAMILY_NAMES <- c(NONE_TYPE_NAME, FLAT_FAMILY_NAME, NORMAL_FAMILY_NAME);
  env$SUPPORTED_COMMON_SCALE_FAMILY_NAMES <- c(NONE_TYPE_NAME, FLAT_FAMILY_NAME, POINT_FAMILY_NAME,
                                               INVGAMMA_FAMILY_NAME, GAMMA_FAMILY_NAME);
  env$UNIVARIATE_FAMILY_NAMES <- c(GAMMA_FAMILY_NAME, INVGAMMA_FAMILY_NAME);
  env$MULTIVARIATE_FAMILY_NAMES <- c(WISHART_FAMILY_NAME, INVWISHART_FAMILY_NAME);
  env$ALL_FAMILY_NAMES <- c(UNIVARIATE_FAMILY_NAMES, MULTIVARIATE_FAMILY_NAMES);
  
  # hyperparameters passed by the user should use these names.
  env$SHAPE_HYPERPARAMETER_NAME <- "shape"
  env$RATE_HYPERPARAMETER_NAME  <- "rate"
  env$SCALE_HYPERPARAMETER_NAME <- "scale"
  env$DEGREES_OF_FREEDOM_HYPERPARAMETER_NAME <- "df"
  env$INVERSE_SCALE_HYPERPARAMETER_NAME <- "inverse.scale"
  env$SD_HYPERPARAMETER_NAME <- "sd";
  env$COVARIANCE_HYPERPARAMETER_NAME <- "cov";
  env$VALUE_HYPERPARAMETER_NAME <- "value";

  # scale-ish things for variance components; 
  env$POSTERIOR_SCALE_OPTION_NAME  <- "posterior.scale"; # compute posterior of sd or cov
  env$DATA_SCALE_OPTION_NAME       <- "data.scale";      # rescale by data
  env$COMMON_SCALE_OPTION_NAME     <- "common.scale";    # variance multiplied by common scale or not
  
  # deprecated
  env$COVARIANCE_SCALE_OPTION_NAME <- "cov.scale";       # is variance to be multiplied by the common scale factor?
  

  
  # useful regular expressions
  # peels the head off of a non-nested comma separated list
  env$commaSeparatedListRegularExpression <-
    "^([^,(]+(?:\\([^)]*\\))?),(.*)$";

  # defaults are specified by not having a twidle
  env$defaultSpecificationRegularExpression <-
    "^\\s*[^~]+\\s*$"

  env$covariancePriorSpecificationPattern <-
    paste("^\\s*(?:([^\\s~]+)\\s*~\\s*)?\\s*(",
          flattenStrings(SUPPORTED_COVARIANCE_TYPE_NAMES, "|"),
          ")(?:\\((.*)\\))?\\s*$", sep="");

  env$unmodeledCoefficientPriorSpecificationPattern <-
    paste("^\\s*(", flattenStrings(SUPPORTED_UNMODELED_COEFFICIENT_FAMILY_NAMES, "|"),
          ")\\s*(?:\\((.*)\\))?\\s*$", sep="");

  env$commonScalePriorSpecificationPattern <-
    paste("^\\s*(", flattenStrings(SUPPORTED_COMMON_SCALE_FAMILY_NAMES, "|"),
          ")\\s*(?:\\((.*)\\))?\\s*$", sep="");
}
