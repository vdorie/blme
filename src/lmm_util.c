#include "lmm.h"

#include "__lmmMerCache.h"

#include "lmer.h"
#include "Syms.h"

#include "matrix.h"

void computeDowndatedDenseCrossproduct(SEXP regression, MERCache* cache, double* target)
{
  const int* dims = DIMS_SLOT(regression);
  
  int numObservations   = dims[n_POS];
  int numModeledCoefs   = dims[q_POS];
  int numUnmodeledCoefs = dims[p_POS];
  
  double* denseDesignMatrix = X_SLOT(regression);
  double* offDiagonalBlockRightFactorization = RZX_SLOT(regression);
  
  int modelIncludesWeights = SXWT_SLOT(regression) != NULL;
  if (modelIncludesWeights) denseDesignMatrix = cache->weightedDenseDesignMatrix;
  
  // downdate X'X
  singleMatrixCrossproduct(denseDesignMatrix, numObservations, numUnmodeledCoefs,
                           target, FALSE, TRIANGLE_TYPE_UPPER);
  
  singleMatrixCrossproductWithUpdate(offDiagonalBlockRightFactorization, numModeledCoefs, numUnmodeledCoefs, -1.0,
                                     target, FALSE, TRIANGLE_TYPE_UPPER);
}
