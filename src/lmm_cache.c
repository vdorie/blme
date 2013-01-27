#include "__lmmMerCache.h"

#include "lmer.h"
#include "blmer_types.h"
#include "Syms.h"

#include "lmm_priors.h"

static void weightDenseComponents(SEXP regression, MERCache* cache);
static commonScaleOptimization_t getCommonScaleOptimizationType(SEXP regression);

void printLMMCache(const MERCache* cache);

MERCache* createLMMCache(SEXP regression)
{
  MERCache* result = (MERCache*) malloc(sizeof(MERCache));
  
  const int* dims = DIMS_SLOT(regression);
  
  int numObservations   = dims[n_POS];
  int numUnmodeledCoefs = dims[p_POS];
  int numModeledCoefs   = dims[q_POS];
  
  // allocation
  result->weightedDenseDesignMatrix      = (double*) malloc(numObservations * numUnmodeledCoefs * sizeof(double));
  result->weightedResponse               = (double*) malloc(numObservations * sizeof(double));
  
  result->downdatedDenseResponseRotation = (double*) malloc(numUnmodeledCoefs * sizeof(double));
  result->downdatedDenseCrossproduct     = (double*) malloc(numUnmodeledCoefs * numUnmodeledCoefs * sizeof(double));
  
  result->unmodeledCoefProjection = (double*) malloc(numUnmodeledCoefs * sizeof(double));
  result->modeledCoefProjection   = (double*) malloc(numModeledCoefs * sizeof(double));
  
  if (result->weightedDenseDesignMatrix == NULL || result->weightedResponse == NULL || result->downdatedDenseResponseRotation == NULL ||
      result->downdatedDenseCrossproduct == NULL || result->unmodeledCoefProjection == NULL || result->modeledCoefProjection == NULL)
    error("Unable to allocate cache, out of memory\n");
  
  weightDenseComponents(regression, result); // also fills in sum( weighted response^2 );
  
  // constants derived from priors
  // flow control
  result->commonScaleOptimization = getCommonScaleOptimizationType(regression);
  
  result->priorDevianceConstantPart = getPriorDevianceConstantPart(regression);
  
  setPriorConstantContributionsToCommonScale(regression, result);
  
  result->mOneExponentialTermVaryingPart = R_NaN;
  result->mTwoExponentialTermVaryingPart = R_NaN;
  result->oneExponentialTermVaryingPart  = R_NaN;
  result->twoExponentialTermVaryingPart  = R_NaN;
  
  result->objectiveFunctionValue = R_NaN;
  result->priorContribution      = R_NaN;
  
  // printLMMCache(result);
  
  return(result);
}

void deleteLMMCache(MERCache* cache)
{
  free(cache->weightedDenseDesignMatrix);
  free(cache->weightedResponse);
  
  free(cache->downdatedDenseResponseRotation);
  free(cache->downdatedDenseCrossproduct);
  free(cache->unmodeledCoefProjection);
  free(cache->modeledCoefProjection);
  free(cache);
}

static void weightDenseComponents(SEXP regression, MERCache* cache)
{
  const int* dims = DIMS_SLOT(regression);
  
  int numObservations   = dims[n_POS];
  int numUnmodeledCoefs = dims[p_POS];
  
  double* sqrtObservationWeight = SXWT_SLOT(regression);
  double* offsets               = OFFSET_SLOT(regression);
  
  const double* denseDesignMatrix   = X_SLOT(regression);
  double* weightedDenseDesignMatrix = cache->weightedDenseDesignMatrix;
  
  const double* response   = Y_SLOT(regression);
  double* weightedResponse = cache->weightedResponse;
  
  cache->responseSumOfSquares = 0.0;
  for (int row = 0; row < numObservations; ++row) {
    double rowWeight = (sqrtObservationWeight ? sqrtObservationWeight[row] : 1.0);
    for (int col = 0; col < numUnmodeledCoefs; ++col) {
      int matrixIndex = row + col * numObservations;
      weightedDenseDesignMatrix[matrixIndex] = denseDesignMatrix[matrixIndex] * rowWeight;
    }
    
    weightedResponse[row] = (response[row] - (offsets ? offsets[row] : 0.0)) * rowWeight;
    cache->responseSumOfSquares += weightedResponse[row] * weightedResponse[row];
  }
}

// the common-scale can be profiled out (easily) under certain circumstances
//
// in general, the objective function looks like:
//   (sigma^2)^(-df/2) * exp(-0.5 * (1 / sigma^2) * stuff)
//
// what determines whether or not it can be profiled is the functional form of the
// exponentiated term. For this, we have the following:
//
// powers: -1 -2  1  2   - estimating equation
// prsnt :  0  X  0  0   - linear in sigma^2 (default scenario)
//          0  X  0  X   - quadratic in sigma^2
//          X  X  0  0   - quadratic in sigma
//          0  X  X  0   - cubic in sigma
//
// everything else is even worse
//
// finally, we also have two trumps. if the common scale has a point prior, that is that.
// In addition, if the unmodeled coefficients aren't placed on the common scale, no polynomial
// equation is possible
static commonScaleOptimization_t getCommonScaleOptimizationType(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  if (!isLinearModel) return(CSOT_NA); // question doesn't apply if is !lmm
  
  int* dims = DIMS_SLOT(regression);
  
  SEXP csPrior     = GET_SLOT(regression, blme_commonScalePriorSym);
  SEXP ucPrior     = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  SEXP cvPriorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  priorType_t csPriorType = PRIOR_TYPE_SLOT(csPrior);
  priorType_t ucPriorType = PRIOR_TYPE_SLOT(ucPrior);
  priorType_t cvPriorType;
  
  priorFamily_t family;
  priorPosteriorScale_t posteriorScale;
  priorCommonScale_t    onCommonScale;
  
  // handle the two trumps first
  if (csPriorType == PRIOR_TYPE_DIRECT &&
      PRIOR_FAMILIES_SLOT(csPrior)[0] == PRIOR_FAMILY_POINT) {
    return(CSOT_NA);
  }
  if (ucPriorType == PRIOR_TYPE_DIRECT) {
    family         = PRIOR_FAMILIES_SLOT(ucPrior)[0];
    onCommonScale  = getCommonScaleBit(PRIOR_SCALES_SLOT(ucPrior)[0]);
    
    if (family != PRIOR_FAMILY_FLAT &&
        (family != PRIOR_FAMILY_GAUSSIAN || !onCommonScale)) return(CSOT_BRUTE_FORCE);
  }
  
  
  // catalog whether or not certain powers are present
  int mOneInExp = 0; // purposefully using 0/1 instead of false/true as we will do some math with them at the end
  int  oneInExp = 0;
  int  twoInExp = 0;
  
  if (csPriorType == PRIOR_TYPE_DIRECT) {
    family         = PRIOR_FAMILIES_SLOT(csPrior)[0];
    posteriorScale = getPosteriorScaleBit(PRIOR_SCALES_SLOT(csPrior)[0]);
    
    if (family == PRIOR_FAMILY_GAMMA) {
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD)
        oneInExp = 1;
      else
        twoInExp = 1;
    } else if (family == PRIOR_FAMILY_INVGAMMA) {
      if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD)
        mOneInExp = 1;
    } else {
      // huh? Shouldn't happen. caught point priors above...
      return(CSOT_BRUTE_FORCE);
    }
  }
  
  // now for the covariance priors; have to loop over factors and over dimensions within
  SEXP stList = GET_SLOT(regression, lme4_STSym);
  int numFactors = dims[nt_POS];
  for (int i = 0; i < numFactors; ++i) {
    SEXP cvPrior  = VECTOR_ELT(cvPriorList, i);
    SEXP stMatrix = VECTOR_ELT(stList, i);
    
    int factorDimension = INTEGER(getAttrib(stMatrix, R_DimSymbol))[0];
    cvPriorType = PRIOR_TYPE_SLOT(cvPrior);
    priorFamily_t* families;
    int* scales;
    
    switch (cvPriorType) {
      case PRIOR_TYPE_DIRECT:
        onCommonScale = getCommonScaleBit(PRIOR_SCALES_SLOT(cvPrior)[0]);
        if (onCommonScale) continue;
        
        family         = PRIOR_FAMILIES_SLOT(cvPrior)[0];
        posteriorScale = getPosteriorScaleBit(PRIOR_SCALES_SLOT(cvPrior)[0]);
        
        switch (family) {
          case PRIOR_FAMILY_GAMMA:
            if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) oneInExp = 1;
            else twoInExp = 1;
            break;
          case PRIOR_FAMILY_INVGAMMA:
            if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) mOneInExp = 1;
            break;
          case PRIOR_FAMILY_WISHART:
            twoInExp = 1;
            break;
          default: break;
        }
        
      case PRIOR_TYPE_CORRELATION:
        families = PRIOR_FAMILIES_SLOT(cvPrior);
        scales   = PRIOR_SCALES_SLOT(cvPrior);
        for (int j = 0; j < factorDimension; ++j) {
          onCommonScale = getCommonScaleBit(scales[j]);
          if (onCommonScale) continue;
          
          family = families[j];
          switch (family) {
            case PRIOR_FAMILY_GAMMA:
              oneInExp = 1;
              break;
            case PRIOR_FAMILY_INVGAMMA:
              mOneInExp = 1;
              break;
            default: break;
          }
        }
        family = families[factorDimension];
        onCommonScale = getCommonScaleBit(scales[factorDimension]);
        if (onCommonScale) continue;
        switch (family) {
          case PRIOR_FAMILY_WISHART:
            oneInExp = 1;
            break;
          case PRIOR_FAMILY_INVWISHART:
            mOneInExp = 1;
            break;
          default: break;
        }
        
      case PRIOR_TYPE_SPECTRAL:
        families = PRIOR_FAMILIES_SLOT(cvPrior);
        scales   = PRIOR_SCALES_SLOT(cvPrior);
        for (int j = 0; j < factorDimension; ++j) {
          onCommonScale = getCommonScaleBit(scales[j]);
          if (onCommonScale) continue;
          
          family = families[j];
          posteriorScale = getPosteriorScaleBit(scales[j]);
          switch (family) {
            case PRIOR_FAMILY_GAMMA:
              if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) oneInExp = 1;
              else twoInExp = 1;
              break;
            case PRIOR_FAMILY_INVGAMMA:
              if (posteriorScale == PRIOR_POSTERIOR_SCALE_SD) mOneInExp = 1;
              break;
            default: break;
          }
        }
      default: break;
    } // switch (cvPriorType)
  } // for (int i = 0; i < numFactors; ++i)
          
    
  int numPowers = mOneInExp + oneInExp + twoInExp;
  if (numPowers == 0) return(CSOT_LINEAR);
  if (numPowers >  1) return(CSOT_BRUTE_FORCE);
  
  if (mOneInExp) return(CSOT_QUADRATIC_SIGMA);
  if (twoInExp)  return(CSOT_QUADRATIC_SIGMA_SQ);
  
  return(CSOT_BRUTE_FORCE);
}

void printLMMCache(const MERCache* cache)
{
  Rprintf("cache: obj = %6f, prr = %6f, yss = %6f, tss = %6f, bss = %6f,\n",
          cache->objectiveFunctionValue, cache->priorContribution, cache->responseSumOfSquares, cache->modeledCoefProjectionSumOfSquares,
          cache->unmodeledCoefProjectionSumOfSquares);
  Rprintf("       tls = %6f, prc = %6f, pdf = %6f, m1c = %6f, m2c = %6f,\n",
          cache->totalSumOfSquares, cache->priorDevianceConstantPart, cache->priorCommonScaleDegreesOfFreedom,
          cache->mOneExponentialTermConstantPart, cache->mTwoExponentialTermConstantPart);
  Rprintf("        1c = %6f,  2c = %6f, m1v = %6f, m2v = %6f,  1v = %6f,  2v = %6f\n",
          cache->oneExponentialTermConstantPart, cache->twoExponentialTermConstantPart,
          cache->mOneExponentialTermVaryingPart, cache->mTwoExponentialTermVaryingPart, cache->oneExponentialTermVaryingPart, cache->twoExponentialTermVaryingPart);
}
