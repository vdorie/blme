#include "commonScalePrior.h"

#include <Rmath.h>
#include <float.h> // used for machine epsilon

#include "Syms.h"
#include "blmer_types.h"
#include "__lmmMerCache.h"

#include "dist_gamma.h"

double getCommonScaleDevianceVaryingPart(SEXP prior, double commonScale)
{
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);
  
  if (priorType != PRIOR_TYPE_DIRECT) return(0.0);
  
  priorFamily_t priorFamily = PRIOR_FAMILIES_SLOT(prior)[0];
  double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  switch (priorFamily) {
    case PRIOR_FAMILY_POINT:
      // shouldn't ever be !=, but can't be sure
      if (fabs(*hyperparameters - commonScale) > DBL_EPSILON) return(INFINITY);
      return(0.0);
      break;
    case PRIOR_FAMILY_GAMMA:
    {
      double shape = hyperparameters[0];
      double rate  = hyperparameters[1];
      priorPosteriorScale_t posteriorScale = getPosteriorScaleBit(PRIOR_SCALES_SLOT(prior)[0]);
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) commonScale = sqrt(commonScale);
      
      return(gammaDevianceVaryingPart(commonScale, shape, rate));
    }
      break;
    case PRIOR_FAMILY_INVGAMMA:
    {
      double shape = hyperparameters[0];
      double scale = hyperparameters[1];
      priorPosteriorScale_t posteriorScale = getPosteriorScaleBit(PRIOR_SCALES_SLOT(prior)[0]);
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) commonScale = sqrt(commonScale);
      
      return(inverseGammaDevianceVaryingPart(commonScale, shape, scale));
    }
      break;
    default:
      break;
  }
  return(0.0);
}

double getCommonScaleDevianceConstantPart(SEXP prior)
{
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);
  
  if (priorType != PRIOR_TYPE_DIRECT) return(0.0);
  
  priorFamily_t priorFamily = PRIOR_FAMILIES_SLOT(prior)[0];
  
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  switch(priorFamily) {
    case PRIOR_FAMILY_GAMMA:
    {
      double shape = hyperparameters[0];
      double rate  = hyperparameters[1];
      
      return(gammaDevianceConstantPart(shape, rate));
    }
      break;
    case PRIOR_FAMILY_INVGAMMA:
    {
      double shape = hyperparameters[0];
      double scale = hyperparameters[1];
    
      return(inverseGammaDevianceConstantPart(shape, scale));
    }
      break;
    default:
      break;
  }
  return(0.0);
}

double getCommonScalePriorDegreesOfFreedom(SEXP prior)
{
  if (PRIOR_TYPE_SLOT(prior) != PRIOR_TYPE_DIRECT) return(0.0);
  
  priorFamily_t priorFamily = PRIOR_FAMILIES_SLOT(prior)[0];
  
  if (priorFamily == PRIOR_FAMILY_POINT) return(0.0);
  
  double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  priorPosteriorScale_t posteriorScale = getPosteriorScaleBit(PRIOR_SCALES_SLOT(prior)[0]);
  
  switch (priorFamily) {
    case PRIOR_FAMILY_GAMMA:
    {
      double shape = hyperparameters[0];
      
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR)
        return(-2.0 * (shape - 1.0));
      else
        return(      -(shape - 1.0));
    }
      break;
    case PRIOR_FAMILY_INVGAMMA:
    {
      double shape = hyperparameters[0];
      
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_VAR)
        return( 2.0 * (shape + 1.0));
      else
        return(        shape + 1.0);
    }
      break;
    default:
      break;
  }
  return(0.0);
}

void setCommonScalePriorExponentialConstantPart(SEXP prior, MERCache* cache)
{
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);
  
  if (priorType != PRIOR_TYPE_DIRECT) return;
  
  priorFamily_t priorFamily = PRIOR_FAMILIES_SLOT(prior)[0];
  priorPosteriorScale_t posteriorScale = getPosteriorScaleBit(PRIOR_SCALES_SLOT(prior)[0]);
  
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  switch (priorFamily) {
    case PRIOR_FAMILY_GAMMA:
    {
      double rate = hyperparameters[1];
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD)
        cache->oneExponentialTermConstantPart += rate;
      else
        cache->twoExponentialTermConstantPart += rate;
    }
      break;
    case PRIOR_FAMILY_INVGAMMA:
    {
      double scale = hyperparameters[1];
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD)
        cache->mOneExponentialTermConstantPart += scale;
      else
        cache->mTwoExponentialTermConstantPart += scale;
    }
      break;
    default:
      break;
  }
}
