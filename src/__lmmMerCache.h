#ifndef __LMM_MER_CACHE__
#define __LMM_MER_CACHE__

typedef enum {
  CSOT_NA,
  CSOT_LINEAR,
  CSOT_QUADRATIC_SIGMA,
  CSOT_QUADRATIC_SIGMA_SQ,
  CSOT_BRUTE_FORCE
} commonScaleOptimization_t;

// supposedly opaque class; not to be mucked around with, and only should be imported into lmm.c

struct _MERCache {
  // anything that needs to be shared between lmm and glmm should come first
  // these *NEED* to be in __glmmMerCache.h and __merCache.h
  double priorDevianceConstantPart;
  
  //
  // anything specific to lmm follows from here
  //
  
  double* weightedDenseDesignMatrix;
  double* weightedResponse;
  
  double* downdatedDenseResponseRotation; // X'y - Lzx * theta half projection
  double* downdatedDenseCrossproduct;     // X'X - LzxRzx
  double* unmodeledCoefProjection;        // Rx^-T (X'y - Lzx * theta half projection)
  double* modeledCoefProjection;          // Lz^-1 A y (A being rotated Z')
  
  double objectiveFunctionValue;
  double priorContribution;
  
  double responseSumOfSquares;
  double modeledCoefProjectionSumOfSquares;
  double unmodeledCoefProjectionSumOfSquares;
  double totalSumOfSquares;

  double mOneExponentialTermConstantPart;
  double mTwoExponentialTermConstantPart;
  double oneExponentialTermConstantPart;
  double twoExponentialTermConstantPart;
  
  // exponential terms; thus far, the only varying part comes from the covariance prior.
  // this little bit is assumed when calculating the covariance deviance using the cache
  double mOneExponentialTermVaryingPart;
  double mTwoExponentialTermVaryingPart;
  double oneExponentialTermVaryingPart;
  double twoExponentialTermVaryingPart;
  
  // in (sigma^2)^(-df/2) * exp(-0.5 * (1 / sigma^2) * stuff), is what influences the 'df' term
  double priorCommonScaleDegreesOfFreedom;
  
  // flow control
  commonScaleOptimization_t commonScaleOptimization;
};

#endif
