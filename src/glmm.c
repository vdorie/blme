#include "glmm.h"

#include "lmer.h"
#include "Syms.h"

#include "unmodeledCoefficientPrior.h"
#include "covariancePrior.h"

#include "__glmmMerCache.h"

MERCache *createGLMMCache(SEXP regression)
{
  MERCache* result = (MERCache *) malloc(sizeof(MERCache));

  const int* dims = DIMS_SLOT(regression);
  
  SEXP unmodeledCoefPrior  = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  
  result->priorDevianceConstantPart =
    getUnmodeledCoefficientDevianceConstantPart(unmodeledCoefPrior, dims[p_POS]) +
    getCovarianceDevianceConstantPart(regression);
    
  return(result);
}

void deleteGLMMCache(MERCache *cache)
{
  free(cache);
}
