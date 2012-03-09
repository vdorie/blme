getCommonScaleDefault <- function(regression, family)
{
  result <- list(family = family);

  fillDefaults <- function(x, ...) x;
  if (result$family == INVGAMMA_FAMILY_NAME) {
    fillDefaults <- getCommonScaleInverseGammaDefaults;
  }
  if (result$family == POINT_FAMILY_NAME) {
    fillDefaults <- getCommonScalePointDefaults;
  }

  result <- fillDefaults(result);
  return(result);
}

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

  prior <- parseCommonScalePrior(regression, options, priorType, callingEnvironment);

  prior <- validateCommonScalePrior(regression, prior);
  prior <- adjustCommonScalePriorScales(regression, prior);

  fields <- getCommonScalePriorFields(regression, prior);
  
  return(createPriorObject(DIRECT_TYPE_NAME, fields));
}

parseCommonScalePrior <- function(regression, specification, family, callingEnvironment)
{
  errorPrefix <- "Error applying prior to common scale: ";
  
  option <- trim(specification);
  if (nchar(option) == 0) return(getCommonScaleDefault(regression, family));
  
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
  
  fillDefaults <- function(x, ...) x;
  if (family == POINT_FAMILY_NAME) {
    prior <- list(family = POINT_FAMILY_NAME);

    prior$value <- getPriorOption(VALUE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    fillDefaults <- getCommonScalePointDefaults;
  } else if (family == INVGAMMA_FAMILY_NAME) {
    prior <- list(family = INVGAMMA_FAMILY_NAME);

    prior$shape <- getPriorOption(SHAPE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$scale <- getPriorOption(SCALE_HYPERPARAMETER_NAME, namedOptionsRef, unnamedOptionsRef);
    prior$posteriorScale <- getPriorOption(POSTERIOR_SCALE_OPTION_NAME, namedOptionsRef, unnamedOptionsRef);

    fillDefaults <- getCommonScaleInverseGammaDefaults;
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

    hyperparameters <- prior$value;
  } else if (prior$family == INVGAMMA_FAMILY_NAME) {
    families <- getEnumOrder(familyEnumeration, INVGAMMA_FAMILY_NAME);
    scales   <- getEnumOrder(scaleEnumeration, prior$posteriorScale);

    shape <- prior$shape;
    if (prior$posteriorScale == SD_SCALE_NAME) {
      # posterior on sd scale with scale = 0, divide shape by 2 to
      # get posterior on var scale
      shape <- shape / 2;
    }
    hyperparameters <- c(shape, prior$scale);
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
    prior$posteriorScale <- defaultCommonScalePointPosteriorScale;

  return(prior);
}

getCommonScaleInverseGammaDefaults <- function(prior)
{
  if (is.null(prior$shape)) {
    prior$shape <- defaultCommonScaleInverseGammaShape;
  }
  
  if (is.null(prior$scale)) {
    prior$scale <- defaultCommonScaleInverseGammaScale;
  }

  if (is.null(prior$posteriorScale)) {
    prior$posteriorScale <- defaultCommonScaleInverseGammaPosteriorScale;
  }

  return(prior);
}

commonScalePriorToString <- function(regression)
{
  if (regression@var.prior@type == getEnumOrder(typeEnumeration, NONE_TYPE_NAME))
    return(character(0));
  
  families <- regression@var.prior@families;
  scales <- regression@var.prior@scales;
  hyperparameters <- regression@var.prior@hyperparameters;

  return(buildStringForFamily(families, scales, hyperparameters, TRUE)$string);
}
