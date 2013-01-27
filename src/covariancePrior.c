#include "covariancePrior.h"

#include <Rmath.h>
#include "Syms.h"

#include "lmer.h"
#include "blmer_types.h"
#include "parameters.h"
#include "util.h"

#include "dist_wishart.h"
#include "dist_gamma.h"

#include "__lmmMerCache.h"

// forward declarations
static double getCorrelationDevianceConstantPart(SEXP prior, int factorDimension);
static double getSpectralDevianceConstantPart(SEXP prior, int factorDimension);
static double getDirectDevianceConstantPart(SEXP prior, int factorDimension);

static double getCorrelationDevianceVaryingPart(SEXP prior, double commonScale, const double *parameters, int factorDimension);
static double getSpectralDevianceVaryingPart(SEXP prior, double commonScale, const double *parameters, int factorDimension);
static double getDirectDevianceVaryingPart(SEXP prior, double commonScale, const double *parameters, int factorDimension);

static void setCorrelationExponentialVaryingParts(SEXP prior, MERCache* cache, const double* parameters, int factorDimension);
static void setSpectralExponentialVaryingParts(SEXP prior, MERCache* cache, const double* parameters, int factorDimension);
static void setDirectExponentialVaryingParts(SEXP prior, MERCache* cache, const double* parameters, int factorDimension);

static double getCorrelationDevianceCommonScaleFreeVaryingPart(SEXP prior, const double *parameters, int factorDimension);
static double getSpectralDevianceCommonScaleFreeVaryingPart(SEXP prior, const double *parameters, int factorDimension);
static double getDirectDevianceCommonScaleFreeVaryingPart(SEXP prior, const double *parameters, int factorDimension);

static double getCorrelationCommonScaleDegreesOfFreedom(SEXP prior, int factorDimension);
static double getSpectralCommonScaleDegreesOfFreedom(SEXP prior, int factorDimension);
static double getDirectCommonScaleDegreesOfFreedom(SEXP prior, int factorDimension);

static void setCorrelationConstraints(double *constraints, int dim);
static void setSpectralConstraints(double *constraints, int dim);
static void setSTConstraints(double *constraints, int dim);

// exported functions
int getNumCovarianceParametersForPrior(priorType_t priorType, int dim)
{
  switch (priorType) {
    case PRIOR_TYPE_CORRELATION:
      return (dim * (dim + 3) / 2);
      break;
    case PRIOR_TYPE_SPECTRAL:
      return (dim * (dim + 1) / 2);
      break;
    case PRIOR_TYPE_DIRECT:
      return (dim * (dim + 1) / 2);
      break;
    case PRIOR_TYPE_NONE:
      return (dim * (dim + 1) / 2);
    default:
      return (0);
      break;
  }
}

double getCovarianceDevianceConstantPart(SEXP regression)
{
  double result = 0.0;
  int numFactors = DIMS_SLOT(regression)[nt_POS];
  
  SEXP stList    = GET_SLOT(regression, lme4_STSym);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  for (int i = 0; i < numFactors; ++i) {
    SEXP st    = VECTOR_ELT(stList, i);
    SEXP prior = VECTOR_ELT(priorList, i);
    
    priorType_t priorType = PRIOR_TYPE_SLOT(prior);
    int factorDimension = INTEGER(getAttrib(st, R_DimSymbol))[0];
    
    switch (priorType) {
      case PRIOR_TYPE_CORRELATION:
        result += getCorrelationDevianceConstantPart(prior, factorDimension);
        break;
      case PRIOR_TYPE_SPECTRAL:
        result += getSpectralDevianceConstantPart(prior, factorDimension);
        break;
      case PRIOR_TYPE_DIRECT:
        result += getDirectDevianceConstantPart(prior, factorDimension);
        break;
      default:
        break;
    }
  }
  
  return(result);
}

double getCovarianceDevianceVaryingPart(SEXP regression, const double* parameters, int* numParametersUsed)
{
  double result = 0.0;
  *numParametersUsed = 0;
  int numParametersForPrior;
  
  const int* dims = DIMS_SLOT(regression);
  
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  double commonScale = 1.0;
  if (isLinearModel) {
    commonScale = DEV_SLOT(regression)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
    commonScale *= commonScale;
  }
  
  int numFactors = dims[nt_POS];
  
  SEXP stList    = GET_SLOT(regression, lme4_STSym);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  for (int i = 0; i < numFactors; ++i) {
    SEXP st    = VECTOR_ELT(stList, i);
    SEXP prior = VECTOR_ELT(priorList, i);
    
    priorType_t priorType = PRIOR_TYPE_SLOT(prior);
    int factorDimension = INTEGER(getAttrib(st, R_DimSymbol))[0];
    
    switch (priorType) {
      case PRIOR_TYPE_CORRELATION:
        result += getCorrelationDevianceVaryingPart(prior, commonScale, parameters, factorDimension);
        break;
      case PRIOR_TYPE_SPECTRAL:
        result += getSpectralDevianceVaryingPart(prior, commonScale, parameters, factorDimension);
        break;
      case PRIOR_TYPE_DIRECT:
        result += getDirectDevianceVaryingPart(prior, commonScale, parameters, factorDimension);
        break;
      default: break;
    }
    
    numParametersForPrior  = getNumCovarianceParametersForPrior(priorType, factorDimension);
    parameters            += numParametersForPrior;
    *numParametersUsed    += numParametersForPrior;
  }

  return(result);
}

void setCovariancePriorCommonScaleExponentialVaryingParts(SEXP regression, MERCache* cache, const double* parameters)
{
  int numFactors = DIMS_SLOT(regression)[nt_POS];
  
  SEXP stList    = GET_SLOT(regression, lme4_STSym);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  for (int i = 0; i < numFactors; ++i) {
    SEXP st    = VECTOR_ELT(stList, i);
    SEXP prior = VECTOR_ELT(priorList, i);
    
    priorType_t priorType = PRIOR_TYPE_SLOT(prior);
    int factorDimension = INTEGER(getAttrib(st, R_DimSymbol))[0];
    
    switch (priorType) {
      case PRIOR_TYPE_CORRELATION:
        setCorrelationExponentialVaryingParts(prior, cache, parameters, factorDimension);
        break;
      case PRIOR_TYPE_SPECTRAL:
        setSpectralExponentialVaryingParts(prior, cache, parameters, factorDimension);
        break;
      case PRIOR_TYPE_DIRECT:
        setDirectExponentialVaryingParts(prior, cache, parameters, factorDimension);
        break;
      default: break;
    }
    
    parameters += getNumCovarianceParametersForPrior(priorType, factorDimension);
  }
}

double getCovarianceDevianceCommonScaleFreeVaryingPart(SEXP regression, const double* parameters)
{
  double result = 0.0;
  
  int numFactors = DIMS_SLOT(regression)[nt_POS];
  
  SEXP stList    = GET_SLOT(regression, lme4_STSym);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  for (int i = 0; i < numFactors; ++i) {
    SEXP st    = VECTOR_ELT(stList, i);
    SEXP prior = VECTOR_ELT(priorList, i);
    
    priorType_t priorType = PRIOR_TYPE_SLOT(prior);
    int factorDimension = INTEGER(getAttrib(st, R_DimSymbol))[0];
    
    switch (priorType) {
      case PRIOR_TYPE_CORRELATION:
        result += getCorrelationDevianceCommonScaleFreeVaryingPart(prior, parameters, factorDimension);
        break;
      case PRIOR_TYPE_SPECTRAL:
        result += getSpectralDevianceCommonScaleFreeVaryingPart(prior, parameters, factorDimension);
        break;
      case PRIOR_TYPE_DIRECT:
        result += getDirectDevianceCommonScaleFreeVaryingPart(prior, parameters, factorDimension);
        break;
      default: break;
    }
    
    parameters += getNumCovarianceParametersForPrior(priorType, factorDimension);
  }
  
  return(result);
}


// common scale has form (sigma^2)^-(df/2) * exp(-0.5 * sigma^-2 * stuff) in the objective function
// adjustment to power of common scale due to prior
double getCovariancePriorCommonScaleDegreesOfFreedom(SEXP regression)
{
  double result = 0.0;
  
  int numFactors = DIMS_SLOT(regression)[nt_POS];
  SEXP stList    = GET_SLOT(regression, lme4_STSym);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  for (int i = 0; i < numFactors; ++i) {
    SEXP st    = VECTOR_ELT(stList, i);
    SEXP prior = VECTOR_ELT(priorList, i);
    
    priorType_t priorType = PRIOR_TYPE_SLOT(prior);
    int factorDimension = INTEGER(getAttrib(st, R_DimSymbol))[0];
    
    switch (priorType) {
      case PRIOR_TYPE_CORRELATION:
        result += getCorrelationCommonScaleDegreesOfFreedom(prior, factorDimension);
        break;
        
      case PRIOR_TYPE_SPECTRAL:
        result += getSpectralCommonScaleDegreesOfFreedom(prior, factorDimension);
        break;
        
      case PRIOR_TYPE_DIRECT:
        result += getDirectCommonScaleDegreesOfFreedom(prior, factorDimension);
        break;
      default:
        break;
    }
  }
  
  
  return(result);
}


void setCovarianceConstraints(SEXP prior, double* boxConstraints, int factorDimension)
{
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);
  switch (priorType) {
    case PRIOR_TYPE_CORRELATION:
      setCorrelationConstraints(boxConstraints, factorDimension);
      break;
      
    case PRIOR_TYPE_SPECTRAL:
      setSpectralConstraints(boxConstraints, factorDimension);
      break;
      
    case PRIOR_TYPE_DIRECT:
      // direct priors are done using the ST decomp
      // setDirectConstraints(boxConstraints, factorDimension);
      setSTConstraints(boxConstraints, factorDimension);
      break;
      
    default:
      setSTConstraints(boxConstraints, factorDimension);
      break;
  }
}

static double getCorrelationDevianceConstantPart(SEXP prior, int factorDimension)
{
  double deviance = 0.0;
  
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  for (int i = 0; i < factorDimension; ++i) {    
    switch (families[i]) {
      case (PRIOR_FAMILY_GAMMA):
      {
        double shape = *hyperparameters++;
        double rate  = *hyperparameters++; // rate is exp( -a * x )
        
        deviance += gammaDevianceConstantPart(shape, rate);
      }
        break;
      case (PRIOR_FAMILY_INVGAMMA):
      {
        double shape = *hyperparameters++;
        double scale = *hyperparameters++; // "scale" here is exp( -a / x ) in the density
        
        deviance += inverseGammaDevianceConstantPart(shape, scale);
      }
        break;
      default:
        break;
    }
  }
  
  switch (families[factorDimension]) {
    case PRIOR_FAMILY_WISHART:
    {
      double degreesOfFreedom    = *hyperparameters++;
      double logScaleDeterminant = *hyperparameters++;
      hyperparameters += factorDimension * factorDimension; // skip scale
      
      deviance += wishartDevianceConstantPart(factorDimension, degreesOfFreedom, logScaleDeterminant);
    }
      break;
    case PRIOR_FAMILY_INVWISHART:
    {
      double degreesOfFreedom           = *hyperparameters++;
      double logInverseScaleDeterminant = *hyperparameters++;
      hyperparameters += factorDimension * factorDimension; // skip inverse.scale
      
      deviance += inverseWishartDevianceConstantPart(factorDimension, degreesOfFreedom, logInverseScaleDeterminant);
    }
      break;
    default:
      break;
  }
  
  return(deviance);
}

static double getCorrelationDevianceVaryingPart(SEXP prior, double commonScale, const double* parameters, int factorDimension)
{
  double deviance = 0.0;
  
  const int* scales             = PRIOR_SCALES_SLOT(prior);
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  priorCommonScale_t    onCommonScale;
  priorPosteriorScale_t posteriorScale;
  
  commonScale = sqrt(commonScale);
  
  for (int i = 0; i < factorDimension; ++i) {
    onCommonScale  = getCommonScaleBit(scales[i]);
    posteriorScale = getPosteriorScaleBit(scales[i]);
    
    switch (families[i]) {
      case (PRIOR_FAMILY_GAMMA):
        {
          double shape = *hyperparameters++;
          double rate  = *hyperparameters++;
          
          double parameter = *parameters++;
          // somewhat counterintuitive, "on common scale" means the prior is placed on
          // the covariance with the common scale divided out, i.e. we don't concern ourselves
          // with it
          if (!onCommonScale) parameter *= commonScale;
          
          if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) parameter *= parameter;
          
          deviance += gammaDevianceVaryingPart(parameter, shape, rate);
        }
        
        break;
        
      case (PRIOR_FAMILY_INVGAMMA):
        {
          double shape = *hyperparameters++;
          double scale = *hyperparameters++; // "scale" here is exp(-a / x) in the density
          
          double parameter = *parameters++;
          
          if (!onCommonScale) parameter *= commonScale;
          
          if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) parameter *= parameter;
          
          deviance += inverseGammaDevianceVaryingPart(parameter, shape, scale);
        }
        
        break;
      default:
        break;
    }
  }
  
  switch (families[factorDimension]) {
    case PRIOR_FAMILY_WISHART:
      {
        double degreesOfFreedom = hyperparameters[0];
        hyperparameters += 2; // skip the log det of the scale
        const double* scaleInverse = hyperparameters;
        hyperparameters += factorDimension * factorDimension;
        
        deviance += wishartDevianceVaryingPart(parameters, factorDimension, degreesOfFreedom, scaleInverse);
      }
      break;
    case PRIOR_FAMILY_INVWISHART:
      {
        double degreesOfFreedom = hyperparameters[0];
        hyperparameters += 2; // skip the log det of the inv scale
        const double* inverseScale = hyperparameters;
        hyperparameters += factorDimension * factorDimension;
        
        deviance += inverseWishartDevianceVaryingPart(parameters, factorDimension, degreesOfFreedom, inverseScale);
      }
      break;
    default:
      break;
  }
    
  parameters += factorDimension * (factorDimension + 1) / 2; // skip an st-style decomp
  
  return (deviance);
}

// The common scale is a bit of a hairy concept here, as there is no unique decomposition.
// We stick with the idea that the center matrix should be really be a correlation matrix,
// and that at some point we might have a better way of applying a prior to it. With that
// in mind, all of the "scale"ness goes with the scale parameters, so all of the common
// scale gets multiplied in there as well. As the decomposition is Scale * coRrelation * Scale,
// we need only multiply each member of the Scale matrix by the square root of the common
// scale to have it multiply out to the real-deal.
//
// This is then slightly complicated by the fact that we permit priors over the scales to be
// over them as variances or standard deviations. As the "scale" is just the standard deviation
// if the correlation is identity, a prior on the variance scale implies a term over the square
// of the common scale - exactly what one might hope.
static void setCorrelationExponentialVaryingParts(SEXP prior, MERCache* cache, const double* parameters, int factorDimension)
{
  const int* scales             = PRIOR_SCALES_SLOT(prior);
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  priorPosteriorScale_t posteriorScale;
  priorCommonScale_t onCommonScale;
  
  for (int i = 0; i < factorDimension; ++i) {
    onCommonScale  = getCommonScaleBit(scales[i]);
    posteriorScale = getPosteriorScaleBit(scales[i]);
    
    switch (families[i]) {
      case (PRIOR_FAMILY_GAMMA):
      {
        double rate = hyperparameters[1];
        hyperparameters += 2; // skip the shape
        
        double parameter = *parameters++;
        
        if (onCommonScale) continue;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) {
          parameter *= parameter;
          // exp(-rate * sigma^2 * Sigma)
          cache->twoExponentialTermVaryingPart += 0.5 * gammaDevianceVaryingExponentialPart(parameter, rate);
        } else {
          // exp(-rate * sigma * Sigma^0.5)
          cache->oneExponentialTermVaryingPart += 0.5 * gammaDevianceVaryingExponentialPart(parameter, rate);
        }
      } 
        break;
        
      case (PRIOR_FAMILY_INVGAMMA):
      {
        double scale = hyperparameters[1];
        hyperparameters += 2; // skip the shape
        
        double parameter = *parameters++;
        
        if (onCommonScale) continue;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) {
          parameter *= parameter;
          // exp(-scale / (sigma^2 * Sigma))
          cache->mTwoExponentialTermVaryingPart += 0.5 * inverseGammaDevianceVaryingExponentialPart(parameter, scale);
        } else {
          // exp(-scale / (sigma * Sigma^0.5))
          cache->mOneExponentialTermVaryingPart += 0.5 * inverseGammaDevianceVaryingExponentialPart(parameter, scale);
        }
      } 
        break;
      default:
        break;
    }
  }
}

static double getCorrelationDevianceCommonScaleFreeVaryingPart(SEXP prior, const double *parameters, int factorDimension)
{
  double result = 0.0;
  
  const int* scales             = PRIOR_SCALES_SLOT(prior);
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  priorCommonScale_t    onCommonScale;
  priorPosteriorScale_t posteriorScale;
  
  for (int i = 0; i < factorDimension; ++i) {
    onCommonScale  = getCommonScaleBit(scales[i]);
    posteriorScale = getPosteriorScaleBit(scales[i]);
    
    switch (families[i]) {
      case (PRIOR_FAMILY_GAMMA):
      {
        double shape = *hyperparameters++;
        double rate  = *hyperparameters++;
        
        double parameter = *parameters++;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) parameter *= parameter;
        
        if (!onCommonScale) {
          // exponential part already taken care of
          result += gammaDevianceVaryingPolynomialPart(parameter, shape);
        } else {
          result += gammaDevianceVaryingPart(parameter, shape, rate);
        }
      } 
        break;
        
      case (PRIOR_FAMILY_INVGAMMA):
      {
        double shape = *hyperparameters++;
        double scale = *hyperparameters++;
        
        double parameter = *parameters++;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) parameter *= parameter;
        
        if (!onCommonScale) {
          result += inverseGammaDevianceVaryingPolynomialPart(parameter, shape);
        } else {
          result += inverseGammaDevianceVaryingPart(parameter, shape, scale);
        }
      } 
        break;
      default:
        break;
    }
  }
  
  return(result);
}

static double getCorrelationCommonScaleDegreesOfFreedom(SEXP prior, int factorDimension) {
  double result = 0.0;
  
  const int* scales             = PRIOR_SCALES_SLOT(prior);
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  double shape;
  priorCommonScale_t   onCommonScale;
  priorPosteriorScale_t posteriorScale;
  
  for (int i = 0; i < factorDimension; ++i) {
    onCommonScale  = getCommonScaleBit(scales[i]);
    posteriorScale = getPosteriorScaleBit(scales[i]);
    
    switch (families[i]) {
      case PRIOR_FAMILY_GAMMA:
        shape = hyperparameters[0];
        hyperparameters += 2; // skip the shape and the rate/scale

        if (onCommonScale) continue;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) {
          // (sigma^2 * param^2) ^ (alpha - 1.0) => (sigma^2)^(-df/2)
          result -= 2.0 * (shape - 1.0); 
        } else {
          // (sigma * param) ^(alpha - 1.0) => (sigma^2)^(-df/2)
          result -=        shape - 1.0;
        }
        break;
        
      case PRIOR_FAMILY_INVGAMMA:
        shape = hyperparameters[0];
        hyperparameters += 2; // skip the shape and the rate/scale
        
        if (onCommonScale) continue;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) {
          result += 2.0 * (shape + 1.0);
        } else {
          result +=        shape + 1.0;
        }
        break;
        
      default: break;
    }
  }
  return(result);
}


static double getSpectralDevianceConstantPart(SEXP prior, int factorDimension)
{
  double deviance = 0.0;
  
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  for (int i = 0; i < factorDimension; ++i) {
    switch (families[i]) {
      case (PRIOR_FAMILY_GAMMA):
      {
        double shape = *hyperparameters++;
        double rate  = *hyperparameters++;
        
        deviance += gammaDevianceConstantPart(shape, rate);
      }
        break;
      case (PRIOR_FAMILY_INVGAMMA):
      {
        double shape = *hyperparameters++;
        double scale = *hyperparameters++;
        
        deviance += inverseGammaDevianceConstantPart(shape, scale);
      }
        break;
        
      default:
        break;
    }
  }
  
  return (deviance);
}


static double getSpectralDevianceVaryingPart(SEXP prior, double commonScale, const double *parameters, int factorDimension)
{
  double deviance = 0.0;
  
  const int* scales             = PRIOR_SCALES_SLOT(prior);
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  priorCommonScale_t onCommonScale;
  priorPosteriorScale_t posteriorScale;

  for (int i = 0; i < factorDimension; ++i) {
    onCommonScale  = getCommonScaleBit(scales[i]);
    posteriorScale = getPosteriorScaleBit(scales[i]);

    switch (families[i]) {
      case (PRIOR_FAMILY_GAMMA):
      {
        double shape = *hyperparameters++;
        double rate  = *hyperparameters++;
        
        // the raw parameter is an eigen value, which is on the scale of a variance
        double parameter = *parameters++;
        
        if (!onCommonScale) parameter *= commonScale;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) parameter = sqrt(parameter);
        
        deviance += gammaDevianceVaryingPart(parameter, shape, rate);
      }
        break;
        
      case (PRIOR_FAMILY_INVGAMMA):
      {
        double shape = *hyperparameters++;
        double scale = *hyperparameters++;
        
        double parameter = *parameters++;
        
        if (!onCommonScale) parameter *= commonScale;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) parameter = sqrt(parameter);
        
        deviance += inverseGammaDevianceVaryingPart(parameter, shape, scale);
      }
        break;

      default:
        break;
    }
  }
  
  return (deviance);
}

static void setSpectralExponentialVaryingParts(SEXP prior, MERCache* cache, const double* parameters, int factorDimension)
{
  const int* scales             = PRIOR_SCALES_SLOT(prior);
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  priorCommonScale_t onCommonScale;
  priorPosteriorScale_t posteriorScale;
  
  for (int i = 0; i < factorDimension; ++i) {
    onCommonScale  = getCommonScaleBit(scales[i]);
    posteriorScale = getPosteriorScaleBit(scales[i]);
    
    switch (families[i]) {
      case (PRIOR_FAMILY_GAMMA):
      {
        double rate = hyperparameters[1];
        hyperparameters += 2; // skip shape
        
        // the raw parameter is an eigen value, which is on the scale of a variance
        double parameter = *parameters++;
        
        if (onCommonScale) continue;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) {
          parameter = sqrt(parameter);
          
          // exp(-rate * sigma * Sigma^0.5)
          cache->oneExponentialTermVaryingPart += 0.5 * gammaDevianceVaryingExponentialPart(parameter, rate);
        } else {
          // exp(-rate * sigma^2 * Sigma)
          cache->twoExponentialTermVaryingPart += 0.5 * gammaDevianceVaryingExponentialPart(parameter, rate);
        }
      }
        break;
        
      case (PRIOR_FAMILY_INVGAMMA):
      {
        double scale = hyperparameters[1];
        hyperparameters += 2; // skip shape
        
        // same as above
        double parameter = *parameters++;
        
        if (onCommonScale) continue;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) {
          parameter = sqrt(parameter);
          // exp(-scale / sigma * Sigma^-0.5)
          cache->mOneExponentialTermVaryingPart += 0.5 * inverseGammaDevianceVaryingExponentialPart(parameter, scale);
        } else {
          // exp(-scale / sigma^2 * Sigma^-1)
          cache->mTwoExponentialTermVaryingPart += 0.5 * inverseGammaDevianceVaryingExponentialPart(parameter, scale);
        }
      }
        break;
        
      default:
        break;
    }
  }
}

static double getSpectralDevianceCommonScaleFreeVaryingPart(SEXP prior, const double *parameters, int factorDimension)
{
  double deviance = 0.0;
  
  const int* scales             = PRIOR_SCALES_SLOT(prior);
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double *hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  priorCommonScale_t onCommonScale;
  priorPosteriorScale_t posteriorScale;
  
  for (int i = 0; i < factorDimension; ++i) {
    onCommonScale  = getCommonScaleBit(scales[i]);
    posteriorScale = getPosteriorScaleBit(scales[i]);
    
    switch (families[i]) {
      case (PRIOR_FAMILY_GAMMA):
      {
        double shape = *hyperparameters++;
        double rate  = *hyperparameters++;
        
        // the raw parameter is an eigen value, which is on the scale of a variance
        double parameter = *parameters++;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) parameter = sqrt(parameter);
        
        if (!onCommonScale) {
          deviance += gammaDevianceVaryingPolynomialPart(parameter, shape);
        } else {
          deviance += gammaDevianceVaryingPart(parameter, shape, rate);
        }
      }
        break;
        
      case (PRIOR_FAMILY_INVGAMMA):
      {
        double shape = *hyperparameters++;
        double scale = *hyperparameters++;
        
        // same as above
        double parameter = *parameters++;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) parameter = sqrt(parameter);
        
        if (!onCommonScale) {
          deviance += inverseGammaDevianceVaryingPolynomialPart(parameter, shape);
        } else {
          deviance += inverseGammaDevianceVaryingPart(parameter, shape, scale);
        }
      }
        break;
        
      default:
        break;
    }
  }
  
  return (deviance);
}


static double getSpectralCommonScaleDegreesOfFreedom(SEXP prior, int factorDimension)
{
  double result = 0.0;
  
  const int* scales             = PRIOR_SCALES_SLOT(prior);
  const priorFamily_t* families = PRIOR_FAMILIES_SLOT(prior);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  double shape;
  priorCommonScale_t onCommonScale;
  priorPosteriorScale_t posteriorScale;
  for (int i = 0; i < factorDimension; ++i) {
    onCommonScale  = getCommonScaleBit(scales[i]);
    posteriorScale = getPosteriorScaleBit(scales[i]);
        
    switch (families[i]) {
      case PRIOR_FAMILY_GAMMA:
        shape = hyperparameters[0];
        hyperparameters += 2; // skip the shape and the rate/scale

        if (onCommonScale) continue;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) {
          result -=        shape - 1.0;
        } else {
          result -= 2.0 * (shape - 1.0);
        }
        break;
        
      case PRIOR_FAMILY_INVGAMMA:
        shape = hyperparameters[0];
        hyperparameters += 2; // skip the shape and the rate/scale

        if (onCommonScale) continue;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) {
          result +=        shape + 1.0;
        } else {
          result += 2.0 * (shape + 1.0);
        }
        break;
      default: break;
    }
  }
  return(result);
}

static double getDirectDevianceConstantPart(SEXP prior, int factorDimension)
{
  double deviance = 0.0;
  
  priorFamily_t family           = PRIOR_FAMILIES_SLOT(prior)[0];
  const double* hyperparameters  = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  switch (family) {
    case (PRIOR_FAMILY_GAMMA):
    {
      double shape = *hyperparameters++;
      double rate  = *hyperparameters++;
      
      deviance += gammaDevianceConstantPart(shape, rate);
    }
      
      break;
    case (PRIOR_FAMILY_INVGAMMA):
    {
      double shape = *hyperparameters++;
      double scale = *hyperparameters++;
      
      deviance += inverseGammaDevianceConstantPart(shape, scale);
    }
      break;
      
    case PRIOR_FAMILY_WISHART:
    {
      double degreesOfFreedom    = *hyperparameters++;
      double logScaleDeterminant = *hyperparameters++;
      // skip the scale inverse
      hyperparameters += factorDimension * factorDimension;
      
      deviance += wishartDevianceConstantPart(factorDimension, degreesOfFreedom, logScaleDeterminant);
    }
      break;
      
    case PRIOR_FAMILY_INVWISHART:
    {
      double degreesOfFreedom           = *hyperparameters++;
      double logInverseScaleDeterminant = *hyperparameters++;
      // skip the inverse scale
      hyperparameters += factorDimension * factorDimension;
      
      deviance += inverseWishartDevianceConstantPart(factorDimension, degreesOfFreedom, logInverseScaleDeterminant);
    }
      break;
      
    default:
      break;
  }
  
  return (deviance);
}


static double getDirectDevianceVaryingPart(SEXP prior, double commonScale, const double* parameters, int factorDimension)
{
  double deviance = 0.0;
  
  int priorScale                = PRIOR_SCALES_SLOT(prior)[0];
  priorFamily_t family          = PRIOR_FAMILIES_SLOT(prior)[0];
  const double *hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  priorCommonScale_t    onCommonScale  = getCommonScaleBit(priorScale);
  priorPosteriorScale_t posteriorScale = getPosteriorScaleBit(priorScale);
  
  double sqrtCommonScale = sqrt(commonScale);
  
  switch (family) {
    case (PRIOR_FAMILY_GAMMA):
      {
        double shape = *hyperparameters++;
        double rate  = *hyperparameters++;
        
        double parameter = parameters[0];
        
        if (!onCommonScale) parameter *= sqrtCommonScale;
        
        // the raw parameter is a standard deviation
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) parameter *= parameter;
        
        deviance += gammaDevianceVaryingPart(parameter, shape, rate);
      }
      
      break;
    case (PRIOR_FAMILY_INVGAMMA):
      {
        double shape = *hyperparameters++;
        double scale = *hyperparameters++;
          
        double parameter = parameters[0];
        
        if (!onCommonScale) parameter *= sqrtCommonScale;
        
        if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) parameter *= parameter;
        
        deviance += inverseGammaDevianceVaryingPart(parameter, shape, scale);
      }
      break;
      
    case PRIOR_FAMILY_WISHART:
      {
        double degreesOfFreedom = hyperparameters[0];
        hyperparameters += 2; // skip the log det of the scale
        const double* scaleInverse = hyperparameters;
        hyperparameters += factorDimension * factorDimension;
        
        if (!onCommonScale) {
          deviance += wishartDevianceVaryingPartWithScale(parameters, factorDimension, degreesOfFreedom, scaleInverse, commonScale);
        } else {
          deviance += wishartDevianceVaryingPart(parameters, factorDimension, degreesOfFreedom, scaleInverse);
        }
      }
      break;
      
    case PRIOR_FAMILY_INVWISHART:
      {
        double degreesOfFreedom = hyperparameters[0];
        hyperparameters += 2; // skip the log det of the inverse scale
        const double* inverseScale = hyperparameters;
        hyperparameters += factorDimension * factorDimension;
        
        if (!onCommonScale) {
          deviance += inverseWishartDevianceVaryingPartWithScale(parameters, factorDimension, degreesOfFreedom, inverseScale, commonScale);
        } else {
          deviance += inverseWishartDevianceVaryingPart(parameters, factorDimension, degreesOfFreedom, inverseScale);
        }
      }
      break;
      
    default:
      break;
  }
  
  return (deviance);
}

static void setDirectExponentialVaryingParts(SEXP prior, MERCache* cache, const double* parameters, int factorDimension)
{
  int priorScale                = PRIOR_SCALES_SLOT(prior)[0];
  priorFamily_t family          = PRIOR_FAMILIES_SLOT(prior)[0];
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  priorCommonScale_t    onCommonScale  = getCommonScaleBit(priorScale);
  priorPosteriorScale_t posteriorScale = getPosteriorScaleBit(priorScale);
  
  if (onCommonScale) return;
  
  switch (family) {
    case (PRIOR_FAMILY_GAMMA):
    {
      double rate = hyperparameters[1];
      hyperparameters += 2; // skip the shape
      
      double parameter = parameters[0];
      
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) {
        parameter *= parameter;
        // exp(-rate * sigma^2 * Sigma)
        cache->twoExponentialTermVaryingPart += 0.5 * gammaDevianceVaryingExponentialPart(parameter, rate);
      } else {
        // exp(-rate * sigma * Sigma^0.5)
        cache->oneExponentialTermVaryingPart += 0.5 * gammaDevianceVaryingExponentialPart(parameter, rate);
      }
    }
      
      break;
    case (PRIOR_FAMILY_INVGAMMA):
    {
      double scale = hyperparameters[1];
      hyperparameters += 2; // skip the shape
      
      double parameter = parameters[0];
      
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) {
        parameter *= parameter;
        cache->mTwoExponentialTermVaryingPart += 0.5 * inverseGammaDevianceVaryingExponentialPart(parameter, scale);
      } else {
        cache->mOneExponentialTermVaryingPart += 0.5 * inverseGammaDevianceVaryingExponentialPart(parameter, scale);
      }
    }
      break;
      
    case PRIOR_FAMILY_WISHART:
    {
      hyperparameters += 2; // skip df and log det scale
      const double* scaleInverse = hyperparameters;
      
      // deviance is -2.0 * log, want just plain -log for "exponential" part
      cache->twoExponentialTermVaryingPart += 0.5 * wishartDevianceVaryingExponentialPart(parameters, factorDimension, scaleInverse);
    }
      break;
      
    case PRIOR_FAMILY_INVWISHART:
    {
      hyperparameters += 2; // skip df and log det inv scale
      const double* inverseScale = hyperparameters;
      
      cache->mTwoExponentialTermVaryingPart += 0.5 * inverseWishartDevianceVaryingExponentialPart(parameters, factorDimension, inverseScale);
    }
      break;
      
    default:
      break;
  }
}

static double getDirectDevianceCommonScaleFreeVaryingPart(SEXP prior, const double *parameters, int factorDimension)
{
  double deviance = 0.0;
  
  int priorScale                = PRIOR_SCALES_SLOT(prior)[0];
  priorFamily_t family          = PRIOR_FAMILIES_SLOT(prior)[0];
  const double *hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  priorCommonScale_t    onCommonScale  = getCommonScaleBit(priorScale);
  priorPosteriorScale_t posteriorScale = getPosteriorScaleBit(priorScale);
  
  switch (family) {
    case (PRIOR_FAMILY_GAMMA):
    {
      double shape = *hyperparameters++;
      double rate  = *hyperparameters++;
      
      double parameter = *parameters++;
      
      // the raw parameter is a standard deviation
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) parameter *= parameter;
      
      if (!onCommonScale) {
        deviance += gammaDevianceVaryingPolynomialPart(parameter, shape);
      } else {
        deviance += gammaDevianceVaryingPart(parameter, shape, rate);
      }
    }
      
      break;
    case (PRIOR_FAMILY_INVGAMMA):
    {
      double shape = *hyperparameters++;
      double scale = *hyperparameters++;
      
      double parameter = *parameters++;
      
      if (parameter == 0.0) return (R_PosInf); // -2 * log(0) = inf
      
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR) parameter *= parameter;
      
      if (!onCommonScale) {
        deviance += inverseGammaDevianceVaryingPolynomialPart(parameter, shape);
      } else {
        deviance += inverseGammaDevianceVaryingPart(parameter, shape, scale);
      }
    }
      break;
      
    case PRIOR_FAMILY_WISHART:
    {
      double degreesOfFreedom = hyperparameters[0];
      hyperparameters += 2; // skip log det of scale
      const double* scaleInverse = hyperparameters;
      
      hyperparameters += factorDimension * factorDimension;
      
      if (!onCommonScale) {
        deviance += wishartDevianceVaryingDeterminantPart(parameters, factorDimension, degreesOfFreedom);
      } else {
        deviance += wishartDevianceVaryingPart(parameters, factorDimension, degreesOfFreedom, scaleInverse);
      }
    }
      break;
      
    case PRIOR_FAMILY_INVWISHART:
    {
      double degreesOfFreedom = hyperparameters[0];
      hyperparameters += 2; // skip log det of inverse scale
      const double* inverseScale = hyperparameters;
      
      hyperparameters += factorDimension * factorDimension;
      
      if (!onCommonScale) {
        deviance += inverseWishartDevianceVaryingDeterminantPart(parameters, factorDimension, degreesOfFreedom);
      } else {
        deviance += inverseWishartDevianceVaryingPart(parameters, factorDimension, degreesOfFreedom, inverseScale);
      }
    }
      break;
      
    default:
      break;
  }
  
  return (deviance);
}

static double getDirectCommonScaleDegreesOfFreedom(SEXP prior, int factorDimension)
{
  int scale                     = PRIOR_SCALES_SLOT(prior)[0];
  priorFamily_t family          = PRIOR_FAMILIES_SLOT(prior)[0];
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  if (getCommonScaleBit(scale) == PRIOR_COMMON_SCALE_TRUE) return(0.0);
  
  switch (family) {
    case PRIOR_FAMILY_GAMMA:
    {
      double shape = hyperparameters[0];
      
      if (getPosteriorScaleBit(scale) == PRIOR_POSTERIOR_SCALE_VAR) {
        return(-2.0 * (shape - 1.0));
      } else {
        return(      -(shape - 1.0));
      }
    } 
      break;
      
    case PRIOR_FAMILY_INVGAMMA:
    {
      double shape = hyperparameters[0];
      
      if (getPosteriorScaleBit(scale) == PRIOR_POSTERIOR_SCALE_VAR) {
        return( 2.0 * (shape + 1.0));
      } else {
        return(       (shape + 1.0));
      }
    }
      break;
      
    case PRIOR_FAMILY_WISHART:
    {
      double degreesOfFreedom = *hyperparameters;
      double d_dim = (double) factorDimension;
      
      return(-d_dim * (degreesOfFreedom - d_dim - 1.0));
    }
      break;
      
    case PRIOR_FAMILY_INVWISHART:
    {
      double degreesOfFreedom = *hyperparameters;
      double d_dim = (double) factorDimension;
      
      return( d_dim * (degreesOfFreedom + d_dim + 1.0));
    }
      break;
      
    default:
      break;
  }
  return(0.0);
}

// constraints
static void setCorrelationConstraints(double* constraints, int dim)
{
  // correlation is scales, diagonal, lower triangle.
  // in theory, we should keep the middle matrix positive definite,
  // but we let our prior handle that
  
  // scales and diagonals
  int numParameters = 2 * dim;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = 0.0;
    *constraints++ = R_PosInf;
  }
  
  // lower triangle of positive definite matrix
  numParameters = dim * (dim - 1) / 2;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = R_NegInf;
    *constraints++ = R_PosInf;
  }
}


static void setSpectralConstraints(double* constraints, int dim)
{
  // spectrals have eigenvalues then angles. The angles should be (-pi, pi],
  // but it is doubtful that the optimization algorithm works on the circle so
  // we extend it to the real line with the understanding that any mode is good
  // enough.
  
  // eigenvalues
  int numParameters = dim;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = 0.0;
    *constraints++ = R_PosInf;
  }
  
  // cayley elements
  numParameters = dim * (dim - 1) / 2;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = R_NegInf;
    *constraints++ = R_PosInf;
  }
}

static void setSTConstraints(double* constraints, int dim)
{
  // default is on the cholesky scale, S first then T
  // S is positive, T is unconstrained
  
  // S
  int numParameters = dim;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = 0.0;
    *constraints++ = R_PosInf;
  }
  
  // T
  numParameters = dim * (dim - 1) / 2;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = R_NegInf;
    *constraints++ = R_PosInf;
  }
}
