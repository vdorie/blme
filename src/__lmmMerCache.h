#ifndef __LMM_MER_CACHE__
#define __LMM_MER_CACHE__

// supposedly opaque class; not to be mucked around with

struct _MERCache {
  double *weightedDenseDesignMatrix;
  double *weightedResponse;
  
  double *downdatedDenseResponseRotation; // X'y - Lzx * theta half projection
  double *downdatedDenseCrossproduct;     // X'X - LzxRzx
  double *unmodeledCoefProjection;        // Rx^-T (X'y - Lzx * theta half projection)
  double *modeledCoefProjection;          // Lz^-1 A y (A being rotated Z')
  
  double responseSumOfSquares;
  double modeledCoefSumOfSquares;
  
  double priorDegreesOfFreedom;
  double priorPenalty;
};

#endif
