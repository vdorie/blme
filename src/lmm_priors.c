#include "lmm_priors.h"

#include "__lmmMerCache.h"

#include "Syms.h"
#include "lmer.h"

#include "covariancePrior.h"
#include "commonScalePrior.h"
#include "unmodeledCoefficientPrior.h"

void setPriorConstantContributionsToCommonScale(SEXP regression, MERCache* cache)
{
  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  SEXP commonScalePrior          = GET_SLOT(regression, blme_commonScalePriorSym);
  
  int numUnmodeledCoefs = DIMS_SLOT(regression)[p_POS];
  
  cache->priorCommonScaleDegreesOfFreedom = 
    getUnmodeledCoefficientPriorCommonScaleDegreesOfFreedom(unmodeledCoefficientPrior, numUnmodeledCoefs) +
    getCommonScalePriorDegreesOfFreedom(commonScalePrior) +
    getCovariancePriorCommonScaleDegreesOfFreedom(regression);
  
  cache->mOneExponentialTermConstantPart = 0.0;
  cache->mTwoExponentialTermConstantPart = 0.0;
  cache->oneExponentialTermConstantPart  = 0.0;
  cache->twoExponentialTermConstantPart  = 0.0;
  
  setCommonScalePriorExponentialConstantPart(GET_SLOT(regression, blme_commonScalePriorSym), cache);
}

void setPriorVaryingContributionsToCommonScale(SEXP regression, MERCache* cache, const double* parameters)
{
  cache->mOneExponentialTermVaryingPart = 0.0;
  cache->mTwoExponentialTermVaryingPart = 0.0;
  cache->oneExponentialTermVaryingPart  = 0.0;
  cache->twoExponentialTermVaryingPart  = 0.0;
  
  setCovariancePriorCommonScaleExponentialVaryingParts(regression, cache, parameters);
}

double getPriorDevianceCommonScaleFreeVaryingPart(SEXP regression, const double* parameters)
{
  return(getCovarianceDevianceCommonScaleFreeVaryingPart(regression, parameters));
}

double getPriorDevianceConstantPart(SEXP regression)
{
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  SEXP commonScalePrior   = GET_SLOT(regression, blme_commonScalePriorSym);

  int numUnmodeledCoefs = DIMS_SLOT(regression)[p_POS];
  
  return(getUnmodeledCoefficientDevianceConstantPart(unmodeledCoefPrior, numUnmodeledCoefs) + 
         getCommonScaleDevianceConstantPart(commonScalePrior) +
         getCovarianceDevianceConstantPart(regression));
}
