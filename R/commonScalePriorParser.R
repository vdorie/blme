parseCommonScalePriorSpecification <- function(regression, specification, callingEnvironment)
{
  if (is.null(specification)) return(createFlatPriorObject());
  
 if (!is.character(specification)) {
    stop("Common scale prior specification is not of character type ",
         "and is unsupported.");
  }

  if (!grepl(commonScalePriorSpecificationPattern, specification, perl=TRUE)) {
    stop(paste("\"", specification, "\" does not parse as a prior specification", sep=""));
  }

  parse <- subSplit(commonScalePriorSpecificationPattern, "\\1\1\\2",
                    specification, perl=TRUE);
  priorType <- parse[1];
  options   <- parse[2];
  if (is.na(options)) options <- "";

  if (priorType == NONE_TYPE_NAME || priorType == FLAT_FAMILY_NAME)
    return (createFlatPriorObject());

  if (priorType == POINT_FAMILY_NAME) {
    prior <- parseCommonScalePointPrior(regression, options, callingEnvironment);
  } else {
    stop("Internal error, please contact the package authors: prior type '",
         priorType, "' unsupported, yet parsed successfully.");
  }

  prior <- validateCommonScalePrior(regression, prior);
  prior <- adjustCommonScalePriorScales(regression, prior);

  fields <- getCommonScalePriorFields(regression, prior);
  
  return(createPriorObject(DIRECT_TYPE_NAME, fields));
}

getCommonScaleDefault <- function(regression, family)
{
  result <- list(family = family);

  fillDefaults <- function(x, ...) x;
  #if (result$family == GAMMA_FAMILY_NAME) {
  #  fillDefaults <- getCommonScaleGammaDefaults;
  #}
  if (result$family == POINT_FAMILY_NAME) {
    fillDefaults <- getCommonScalePointDefaults;
  }

  result <- fillDefaults(result);
  return(result);
}

parseCommonScalePointPrior <- function(regression, specification, callingEnvironment)
{
  errorPrefix <- "Error applying prior to common scale: ";
  
  option <- trim(specification);
  if (nchar(option) == 0) return(getCommonScaleDefault(regression, POINT_FAMILY_NAME));
  
  optionsList <- eval(parse(text=paste("list(", option, ")", sep="")),
                      envir=callingEnvironment);
  
  if (is.null(names(optionsList))) {
    namedOptions   <- list();
    unnamedOptions <- optionsList;
  } else {
    namedOptions   <- optionsList[names(optionsList) != ""];
    unnamedOptions <- optionsList[names(optionsList) == ""];
  }

  # references
  namedOptionsRef   <- list(env = sys.frame(sys.nframe()), var = "namedOptions");
  unnamedOptionsRef <- list(env = sys.frame(sys.nframe()), var = "unnamedOptions");

  family <- POINT_FAMILY_NAME; # would otherwise be passed in, but nothing else is supported
  
  fillDefaults <- function(x, ...) x;
  if (family == POINT_FAMILY_NAME) {
    prior <- list(family = POINT_FAMILY_NAME);

    prior$value <- getPriorOption(VALUE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    fillDefaults <- getCommonScalePointDefaults;
  }

  if (length(namedOptions) > 0) {
    warning("Unrecognized prior option(s) for ", family, " family: ",
            toString(names(namedOptions)), ".");
  }
  if (length(unnamedOptions) > 0) {
    warning("Extra option(s) for ", family, " family: ", toString(unnamedOptions), ".");
  }

  prior <- fillDefaults(prior);

  return(prior);
}

getCommonScalePriorFields <- function(regression, prior)
{
  numUnmodeledCoefs <- regression@dims[["p"]];
  
  families <- integer(0);
  scales   <- integer(0);
  hyperparameters <- double(0);
  
  if (prior$family == POINT_FAMILY_NAME) {
    families <- getEnumOrder(familyEnumeration, POINT_FAMILY_NAME);
    scales   <- getEnumOrder(scaleEnumeration, prior$posteriorScale);

    if (length(scales) == 0) {
      stop("Unable to recognize covariance scale '", prior[["covarianceScale"]], "'.");
    }

    hyperparameters <- prior$value;
  }

  return (list(families        = families,
               scales          = scales,
               hyperparameters = hyperparameters));
}

getCommonScalePointDefaults <- function(prior)
{
  if (is.null(prior$value)) {
    prior$value <- defaultCommonScalePointPriorValue;
  }

  if (is.null(prior$posteriorScale))
    prior$posteriorScale <- defaultCommonScalePointPriorScale;

  return(prior);
}

commonScalePriorToString <- function(regression)
{
  if (regression@var.prior@type == getEnumOrder(typeEnumeration, NONE_TYPE_NAME))
    return(character(0));
  
  families <- regression@var.prior@families;
  scales <- regression@var.prior@scales;
  hyperparameters <- regression@var.prior@hyperparameters;

  return(buildStringForFamily(families, scales, hyperparameters));
}
