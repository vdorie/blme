#include "lmm.h"

#include <Matrix.h>       // cholesky matrix types
#include <R_ext/Lapack.h> // native matrix stuff not worth putting in "matrix.h"

#include "Syms.h"
#include "lmer.h"
#include "lmer_common.h"
#include "util.h"
#include "matrix.h"
#include "blmer.h"
#include "unmodeledCoefficientPrior.h"

#include "common_inlines.h"

extern cholmod_common cholmodCommon;

struct _MERCache {
  double *weightedDenseDesignMatrix;
  double *weightedResponse;
  
  double *downdatedDenseResponseRotation; // X'y - Lzx * theta half projection
  double *downdatedDenseCrossproduct;     // X'X - LzxRzx
  double *unmodeledCoefProjection;    // Rx^-T (X'y - Lzx * theta half projection)
  double *modeledCoefProjection;      // Lz^-1 A y (A being rotated Z')
  
  double responseSumOfSquares;
  double modeledCoefSumOfSquares;
};


#define COMMON_SCALE_OPTIMIZATION_TOLERANCE 1.0e-10

// used by updateWeights
//static void updateResidualWeights(SEXP regression);
//static void updateObservationWeights(SEXP regression);
static void weightSparseDesignMatrix(SEXP regression);
static void weightResponseAndDenseDesignMatrix(SEXP regression, MERCache *cache);

// used by updateAugmentedDesignMatrixFactorizations
static void updateUpperLeftCholeskyBlock(SEXP regression);
static void updateOffDiagonalCholeskyBlock(SEXP regression, MERCache *cache);
static void updateLowerRightCholeskyBlock(SEXP regression, MERCache *cache);

// used to update the common scale when it can't be profiled out
static void
getDerivatives(SEXP regression, MERCache *cache, double *firstDerivative, double *secondDerivative);
static void
updateRegressionForNewCommonScale(SEXP regression, MERCache *cache);

MERCache *createLMMCache(SEXP regression)
{
  int *dims = DIMS_SLOT(regression);
  
  int numObservations   = dims[n_POS];
  int numUnmodeledCoefs = dims[p_POS];
  int numModeledCoefs   = dims[q_POS];
  
  MERCache *result = (MERCache *) malloc(sizeof(MERCache));
  result->weightedDenseDesignMatrix      = (double *) malloc(numObservations * numUnmodeledCoefs * sizeof(double));
  result->weightedResponse               = (double *) malloc(numObservations * sizeof(double));
  
  result->downdatedDenseResponseRotation = (double *) malloc(numUnmodeledCoefs * sizeof(double));
  result->downdatedDenseCrossproduct     = (double *) malloc(numUnmodeledCoefs * numUnmodeledCoefs * sizeof(double));
  
  result->unmodeledCoefProjection = (double *) malloc(numUnmodeledCoefs * sizeof(double));
  result->modeledCoefProjection   = (double *) malloc(numModeledCoefs * sizeof(double));
  
  double *sqrtObservationWeight = SXWT_SLOT(regression);
  double *offsets               = OFFSET_SLOT(regression);
    
  const double *denseDesignMatrix   = X_SLOT(regression);
  double *weightedDenseDesignMatrix = result->weightedDenseDesignMatrix;
  
  const double *response   = Y_SLOT(regression);
  double *weightedResponse = result->weightedResponse;
  
  result->responseSumOfSquares = 0.0;
  for (int row = 0; row < numObservations; ++row) {
    double rowWeight = (sqrtObservationWeight ? sqrtObservationWeight[row] : 1.0);
    for (int col = 0; col < numUnmodeledCoefs; ++col) {
      int matrixIndex = row + col * numObservations;
      weightedDenseDesignMatrix[matrixIndex] = denseDesignMatrix[matrixIndex] * rowWeight;
    }
      
    weightedResponse[row] = (response[row] - (offsets ? offsets[row] : 0.0)) * rowWeight;
    result->responseSumOfSquares += weightedResponse[row] * weightedResponse[row];
  }
  
  return(result);
}

void deleteLMMCache(MERCache *cache)
{
  free(cache->weightedDenseDesignMatrix);
  free(cache->weightedResponse);
  
  free(cache->downdatedDenseResponseRotation);
  free(cache->downdatedDenseCrossproduct);
  free(cache->unmodeledCoefProjection);
  free(cache->modeledCoefProjection);
  free(cache);
}

double lmmCalculateDeviance(SEXP regression, MERCache *cache)
{
  int    *dims      = DIMS_SLOT(regression);
  double *deviances = DEV_SLOT(regression);
  
  rotateSparseDesignMatrix(regression);
#ifdef PRINT_TRACE
  CHM_SP A = A_SLOT(regression);
  R_CheckStack();
#endif
  DEBUG_PRINT_ARRAY("A.x", A->x, ((int *) A->p)[dims[n_POS]] > 10 ? 10 : ((int *) A->p)[dims[n_POS]]);
  
  updateWeights(regression, cache);
  updateAugmentedDesignMatrixFactorizations(regression, cache);
  
  if (canProfileCommonScale(regression)) {
    calculateProjections(regression, cache);
    calculatePenalizedWeightedResidualSumOfSquaresFromProjections(regression, cache);
    rotateProjections(regression, cache);
    
    updateDeviance(regression);
    profileCommonScale(regression);
  } else {
    calculateProjections(regression, cache);
    calculatePenalizedWeightedResidualSumOfSquaresFromProjections(regression, cache);
    
    double currCommonScale = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
    double prevCommonScale;
  
    // Rprintf("starting with scale: %f\n", currCommonScale);
    do {
      prevCommonScale = currCommonScale;
      
      // implictly updates projections and pwrss
      currCommonScale = performOneStepOfNewtonsMethodForCommonScale(regression, cache);
      // Rprintf("  updated to scale: %f\n", currCommonScale);
    } while (fabs(currCommonScale - prevCommonScale) >= COMMON_SCALE_OPTIMIZATION_TOLERANCE);
    
    rotateProjections(regression, cache);
    updateDeviance(regression);
  }
  
  return (deviances[dims[isREML_POS] ? REML_POS : ML_POS]);
}

double lmmApproximateDeviance(SEXP regression, MERCache *cache)
{
  int    *dims      = DIMS_SLOT(regression);
  double *deviances = DEV_SLOT(regression);
  
  rotateSparseDesignMatrix(regression);
  updateWeights(regression, cache);
  updateAugmentedDesignMatrixFactorizations(regression, cache);
  
  calculateProjections(regression, cache);
  calculatePenalizedWeightedResidualSumOfSquaresFromProjections(regression, cache);
  rotateProjections(regression, cache);
  
  updateDeviance(regression);
  // note, we do not do: profileCommonScale(regression);

  return (deviances[dims[isREML_POS] ? REML_POS : ML_POS]);
}

void updateWeights(SEXP regression, MERCache *cache)
{
  // see function definitions as to why these are commented out
  // updateResidualWeights(regression);
  // updateObservationWeights(regression);
  
  weightSparseDesignMatrix(regression);
  weightResponseAndDenseDesignMatrix(regression, cache);
}

void updateAugmentedDesignMatrixFactorizations(SEXP regression, MERCache *cache)
{
  updateUpperLeftCholeskyBlock(regression);
  
#ifdef PRINT_TRACE
  int *dims = DIMS_SLOT(regression);
  CHM_FR L = L_SLOT(regression);
  R_CheckStack();
#endif
  
  DEBUG_PRINT_ARRAY("L.x", L->x, ((int *) L->p)[dims[q_POS]] > 10 ? 10 : ((int *) L->p)[dims[q_POS]]);
  
  updateOffDiagonalCholeskyBlock(regression, cache);
  
  DEBUG_PRINT_ARRAY("RZX", RZX_SLOT(regression), dims[p_POS] * dims[q_POS] > 10 ? 10 : dims[p_POS] * dims[q_POS]);
  
  updateLowerRightCholeskyBlock(regression, cache);
  
  DEBUG_PRINT_ARRAY("RX ", RX_SLOT(regression), dims[p_POS] * dims[p_POS] > 10 ? 10 : dims[p_POS] * dims[p_POS]);
}

void updateDeviance(SEXP regression)
{
  int    *dims      = DIMS_SLOT(regression);
  double *deviances = DEV_SLOT(regression);
  
  int numModeledCoefs = dims[q_POS];
  
  deviances[usqr_POS] = getSumOfSquares(U_SLOT(regression), numModeledCoefs);
  // weighted residual sum of squares is the penalized version, minus the penalty
  deviances[wrss_POS] = deviances[disc_POS] = deviances[pwrss_POS] - deviances[usqr_POS];
  
  double   MLDegreesOfFreedom = (double)  dims[n_POS];
  double REMLDegreesOfFreedom = (double) (dims[n_POS] - dims[p_POS]);
  
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  if (PRIOR_TYPE_SLOT(unmodeledCoefPrior) == PRIOR_TYPE_DIRECT &&
      PRIOR_FAMILIES_SLOT(unmodeledCoefPrior)[0] == PRIOR_FAMILY_GAUSSIAN &&
      PRIOR_SCALES_SLOT(unmodeledCoefPrior)[0] == PRIOR_SCALE_COMMON)
  {
      MLDegreesOfFreedom += (double) dims[p_POS];
    REMLDegreesOfFreedom += (double) dims[p_POS];
  }

  // see attached math for how this is the deviance
  if (canProfileCommonScale(regression)) {    
    deviances[ML_POS] = deviances[ldL2_POS] +
        MLDegreesOfFreedom * (1.0 + log(deviances[pwrss_POS]) + log(2.0 * PI /   MLDegreesOfFreedom));
    deviances[REML_POS] = deviances[ldL2_POS] + deviances[ldRX2_POS] +
      REMLDegreesOfFreedom * (1.0 + log(deviances[pwrss_POS]) + log(2.0 * PI / REMLDegreesOfFreedom));
  } else {
    double sigma_sq = deviances[sigmaML_POS] * deviances[sigmaML_POS];
    deviances[ML_POS] = deviances[ldL2_POS] + 
      MLDegreesOfFreedom * log(2.0 * PI * sigma_sq) + deviances[pwrss_POS] / sigma_sq;
    
    sigma_sq = deviances[sigmaREML_POS] * deviances[sigmaREML_POS];
    deviances[REML_POS] = deviances[ldL2_POS] + deviances[ldRX2_POS] +
      REMLDegreesOfFreedom * log(2.0 * PI * sigma_sq) + deviances[pwrss_POS] / sigma_sq;
  }
}

void profileCommonScale(SEXP regression)
{
  int    *dims      = DIMS_SLOT(regression);
  double *deviances = DEV_SLOT(regression);
  
  int numObservations = dims[n_POS];

  double   MLDegreesOfFreedom = (double)  dims[n_POS];
  double REMLDegreesOfFreedom = (double) (dims[n_POS] - dims[p_POS]);
  
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  if (PRIOR_TYPE_SLOT(unmodeledCoefPrior) == PRIOR_TYPE_DIRECT &&
      PRIOR_FAMILIES_SLOT(unmodeledCoefPrior)[0] == PRIOR_FAMILY_GAUSSIAN &&
      PRIOR_SCALES_SLOT(unmodeledCoefPrior)[0] == PRIOR_SCALE_COMMON)
  {
      MLDegreesOfFreedom += (double) dims[p_POS];
    REMLDegreesOfFreedom += (double) dims[p_POS];
  }
  
  double *sqrtResidualWeight = SRWT_SLOT(regression);
  
  deviances[sigmaML_POS]   = sqrt(deviances[pwrss_POS] /
                                  (sqrtResidualWeight ? getSumOfSquares(sqrtResidualWeight, numObservations) : MLDegreesOfFreedom));
  deviances[sigmaREML_POS] = deviances[sigmaML_POS] * sqrt(MLDegreesOfFreedom / REMLDegreesOfFreedom);
}

double performOneStepOfNewtonsMethodForCommonScale(SEXP regression, MERCache *cache)
{
  double *deviances = DEV_SLOT(regression);
  int *dims = DIMS_SLOT(regression);
  
  double firstDerivative;
  double secondDerivative;
  
  getDerivatives(regression, cache, &firstDerivative, &secondDerivative);
  
  double oldCommonScale = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  double commonScaleDelta = firstDerivative / secondDerivative;
  double newCommonScale = oldCommonScale - commonScaleDelta;
  
  if (newCommonScale < 0.0) {
    newCommonScale = oldCommonScale / 2.0; // revert to bissection
  } else if (commonScaleDelta < 0.0 && secondDerivative > 0.0) {
    // There should be a convex region that extends from 0 up to a point.
    // When we're past that point, updating to larger values of
    // the scale parameter is a bad idea (should asymptote somewhere
    // below a root). We can identify this region by the change in
    // sign of the second derivative (that the log-likelihood goes
    // to -Inf at 0 and the convex region gives us that the sign
    // should be positive).
    //
    // When this happens, we need to move inward, not outward. Could
    // perhaps do something scaled to the second derivative, but
    // this seems as reasonable as anything else.
    newCommonScale = oldCommonScale / 2.0;
  }
  
  deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS] = newCommonScale;
  
  updateRegressionForNewCommonScale(regression, cache);
  
  return (deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS]);
}

void calculateProjectionsForSingleArgumentAnova(SEXP regression,
                                                double *modeledCoefProjection,
                                                double *unmodeledCoefProjection)
{
  MERCache *cache = createLMMCache(regression);
  
  // replace where the cache points to save a little bit of time
  // as opposed to copying from the temp
  double *cacheModeledCoefProjection = cache->modeledCoefProjection;
  double *cacheUnmodeledCoefProjection = cache->unmodeledCoefProjection;
  
  cache->modeledCoefProjection   = modeledCoefProjection;
  cache->unmodeledCoefProjection = unmodeledCoefProjection;
  
  updateWeights(regression, cache);
  calculateProjections(regression, cache);
  calculatePenalizedWeightedResidualSumOfSquaresFromProjections(regression, cache);

  cache->modeledCoefProjection   = cacheModeledCoefProjection;
  cache->unmodeledCoefProjection = cacheUnmodeledCoefProjection;
  
  deleteLMMCache(cache);
}

// This function actually should never be called, because it is (should be?)
// entirely superfluous for a linear model.
/* static void updateResidualWeights(SEXP regression)
{
  int    *dims       = DIMS_SLOT(regression);
  double *deviances  = DEV_SLOT(regression);
  
  int numObservations = dims[n_POS];
  
  double *variances    = VAR_SLOT(regression);
  double *priorWeights = PWT_SLOT(regression);
  
  if (variances == NULL && priorWeights == NULL) return;
  
  double *sqrtResidualWeights = SRWT_SLOT(regression);
  
  double *residuals = RESID_SLOT(regression);
  double *response  = Y_SLOT(regression);
  double *mu        = MU_SLOT(regression); // current estimate of E(y)
  
  deviances[wrss_POS] = 0.0;
    
  for (int i = 0; i < numObservations; ++i) {
    double priorWeight = (priorWeights ? priorWeights[i] : 1.0);
    double variance    = (variances    ? variances[i]    : 1.0);
    
    sqrtResidualWeights[i] = sqrt(priorWeight / variance);
    residuals[i]           = sqrtResidualWeights[i] * (response[i] - mu[i]);
    
    deviances[wrss_POS] += residuals[i] * residuals[i];
  }
} */

// as with updating the residual weights, for a linear model if the
// weights exist at all, they are unchanged across every iteration
/* static void updateObservationWeights(SEXP regression)
{
  int *dims = DIMS_SLOT(regression);
  
  int numObservations = dims[n_POS];
  int numReplicationsOfObservations = dims[s_POS]; // definitely needs a better name
  
  double *sqrtResidualWeight = SRWT_SLOT(regression);
  
  // the model row weight is further scaled by the derivative of the mean
  // w.r.t. the linear predictor
  double *sqrtObservationWeight = SXWT_SLOT(regression); // W^1/2 * diag(d mu / d eta)
  double *gradient = V_SLOT(regression);
  double *muEta    = MUETA_SLOT(regression);
  
  if (sqrtObservationWeight != NULL) {			// Update sXwt - matrix of model row weights 
    for (int j = 0; j < numReplicationsOfObservations; ++j) {  // numReps == 1 unless NLMM 
	    for (int i = 0; i < numObservations; ++i) {
        int weightIndex = i + j * numObservations;
        
        sqrtObservationWeight[weightIndex] =
          (sqrtResidualWeight ? sqrtResidualWeight[i] : 1.0) *
		      (muEta              ? muEta[i]              : 1.0) *
          (gradient           ? gradient[weightIndex] : 1.0); // gradient is NULL unless NLMM
	    }
    }
  }
} */

static void weightSparseDesignMatrix(SEXP regression)
{
  double *sqrtObservationWeight = SXWT_SLOT(regression); // W^1/2 * diag(d mu / d eta)
  
  if (sqrtObservationWeight == NULL) return;

  int *dims = DIMS_SLOT(regression);
  
  int numObservations = dims[n_POS];

  CHM_SP rotatedSparseDesignMatrix = A_SLOT(regression);
  R_CheckStack();
  
  // int *nonZeroRowIndices = (int *) rotatedSparseDesignMatrix->i;
  int *indicesForColumn  = (int *) rotatedSparseDesignMatrix->p;
  double *values = (double *) rotatedSparseDesignMatrix->x;
  
  double *cMatrixExplicitStorage = Cx_SLOT(regression); // used when Cm doesn't exist

  // cMatrixExplicitStorage exists only when numReplicationsOfObservations = 1, which should be for all GLMM and LMM
  if (cMatrixExplicitStorage == NULL) error("For an lmm with weights, Cx is required to exist but has not been defined.");
  
  for (int col = 0; col < numObservations; ++col) {
    // since A/C are (rotated) Z, transposed, the columns of A/C correspond to observations
    for (int sparseIndex = indicesForColumn[col]; sparseIndex < indicesForColumn[col + 1]; ++sparseIndex) {
      cMatrixExplicitStorage[sparseIndex] = values[sparseIndex] * sqrtObservationWeight[col];
    }
  }
}

static void weightResponseAndDenseDesignMatrix(SEXP regression, MERCache *cache)
{
  int *dims = DIMS_SLOT(regression);
  
  int numObservations   = dims[n_POS];
  int numUnmodeledCoefs = dims[p_POS];
  
  double *response         = Y_SLOT(regression);
  double *weightedResponse = cache->weightedResponse;
  double *offsets          = OFFSET_SLOT(regression);
  
  for (int i = 0; i < numObservations; ++i) weightedResponse[i] = response[i] - (offsets ? offsets[i] : 0.0);
  
  double *sqrtObservationWeight = SXWT_SLOT(regression);
  if (sqrtObservationWeight == NULL) return;
  
  double *denseDesignMatrix         = X_SLOT(regression);
  double *weightedDenseDesignMatrix = cache->weightedDenseDesignMatrix;
  cache->responseSumOfSquares = 0.0;
  
  for (int i = 0; i < numObservations; ++i) {
    weightedResponse[i] *= sqrtObservationWeight[i];
    cache->responseSumOfSquares += weightedResponse[i] * weightedResponse[i];
    
    for (int j = 0; j < numUnmodeledCoefs; ++j) {
      int elementIndex = i + j * numObservations;
      
      weightedDenseDesignMatrix[elementIndex] = sqrtObservationWeight[i] * denseDesignMatrix[elementIndex];
    }
  }
}

static void updateUpperLeftCholeskyBlock(SEXP regression)
{
  CHM_SP rotatedSparseDesignMatrix       = A_SLOT(regression);
  CHM_FR upperLeftBlockLeftFactorization = L_SLOT(regression);
  R_CheckStack();
  
  int modelIncludesWeights = SXWT_SLOT(regression) != NULL;
  
  // replace A with C if needed, as previously computed
  if (modelIncludesWeights) {
    double *cMatrixExplicitStorage = Cx_SLOT(regression);
    
    if (cMatrixExplicitStorage == NULL) {
      error("For an lmm with weights, Cx is required to exist but has not been defined.");
    }
    rotatedSparseDesignMatrix->x = (void *) cMatrixExplicitStorage;
  }
  rotatedSparseDesignMatrix->stype = 0;

  
  double d_one[] = { 1.0, 0 };
  
  // factorize I + A'A, or I + C'C if it exists, where C = AW^(1/2), and A = ST'Z'
  if (!M_cholmod_factorize_p(rotatedSparseDesignMatrix, d_one, NULL, 0 /*fsize*/,
                             upperLeftBlockLeftFactorization, &cholmodCommon)) {
    error("cholmod_factorize_p failed: status %d, minor %d from ncol %d",
          cholmodCommon.status,
          upperLeftBlockLeftFactorization->minor,
          upperLeftBlockLeftFactorization->n);
  }
  
  double *deviances = DEV_SLOT(regression);
  deviances[ldL2_POS] = M_chm_factor_ldetL2(upperLeftBlockLeftFactorization);
}

static void updateOffDiagonalCholeskyBlock(SEXP regression, MERCache *cache)
{
  int *dims = DIMS_SLOT(regression);
  int numUnmodeledCoefs = dims[p_POS];

  double *denseDesignMatrix = X_SLOT(regression);
  double *offDiagonalBlockRightFactorization = RZX_SLOT(regression);
    
  CHM_SP rotatedSparseDesignMatrix       = A_SLOT(regression);
  CHM_FR upperLeftBlockLeftFactorization = L_SLOT(regression);
  R_CheckStack();

  int modelIncludesWeights = SXWT_SLOT(regression) != NULL;
  
  // replace model matrices by their weighted versions if necessary
  if (modelIncludesWeights) {
    denseDesignMatrix = cache->weightedDenseDesignMatrix;
    
    double *cMatrixExplicitStorage = Cx_SLOT(regression);
    if (cMatrixExplicitStorage == NULL) {
      error("For an lmm with weights, Cx is required to exist but has not been defined.");
    }
    
    rotatedSparseDesignMatrix->x = (void *) cMatrixExplicitStorage;
  }
  
  // RZX (temp) = PAX, or
  // temp = Perm * (weighted) rot.Z' * (weighted) X
  multiplyWithPermutation(offDiagonalBlockRightFactorization,
                          (int *) upperLeftBlockLeftFactorization->Perm,
                          rotatedSparseDesignMatrix,
                          denseDesignMatrix,
                          numUnmodeledCoefs);
                          
  // solve L %*% RZX = PAW^{1/2}GHX, or
  // prodFactor = L^-1 * Perm * rot.Z' * X
  // L being left factor of upper left quadrant of augmented design, rot.Z * rot.Z' + I
  solveSparseCholeskySystem(CHOLMOD_L, upperLeftBlockLeftFactorization,
                            offDiagonalBlockRightFactorization, numUnmodeledCoefs,
                            offDiagonalBlockRightFactorization);
}

static void updateLowerRightCholeskyBlock(SEXP regression, MERCache *cache)
{
  int *dims = DIMS_SLOT(regression);
  
  int numObservations   = dims[n_POS];
  int numModeledCoefs   = dims[q_POS];
  int numUnmodeledCoefs = dims[p_POS];
  
  double *denseDesignMatrix = X_SLOT(regression);
  
  double *offDiagonalBlockRightFactorization = RZX_SLOT(regression);
  double *lowerRightBlockRightFactorization  = RX_SLOT(regression);
  
  
  int modelIncludesWeights = SXWT_SLOT(regression) != NULL;
  
  if (modelIncludesWeights) denseDesignMatrix = cache->weightedDenseDesignMatrix;

  
  /* downdate X'X and factor  */
  // compute X'X
  singleMatrixCrossproduct(denseDesignMatrix, numObservations, numUnmodeledCoefs,
                           lowerRightBlockRightFactorization, FALSE, TRIANGLE_TYPE_UPPER);

  singleMatrixCrossproductWithUpdate(offDiagonalBlockRightFactorization, numModeledCoefs, numUnmodeledCoefs, -1.0,
                                     lowerRightBlockRightFactorization, FALSE, TRIANGLE_TYPE_UPPER);

  if (!canProfileCommonScale(regression)) {
    Memcpy(cache->downdatedDenseCrossproduct, (const double *) lowerRightBlockRightFactorization,
           numUnmodeledCoefs * numUnmodeledCoefs);
  }
  
  addGaussianContributionToDenseBlock(regression, lowerRightBlockRightFactorization);
  
  int choleskyResult = getDenseCholeskyDecomposition(lowerRightBlockRightFactorization, numUnmodeledCoefs,
                                                     TRIANGLE_TYPE_UPPER);

  if (choleskyResult > 0) error("Leading minor %d of downdated X'X is not positive definite.", choleskyResult);
  if (choleskyResult < 0) error("Illegal argument %d to cholesky decomposition (dpotrf).", -choleskyResult);
  
  double *deviances = DEV_SLOT(regression);
  // accumulate log(det(RX)^2)
  // matrix is triangular
  deviances[ldRX2_POS] = 0.0;
  for (int j = 0; j < numUnmodeledCoefs; ++j) {
    deviances[ldRX2_POS] += 2.0 * log(lowerRightBlockRightFactorization[j * (numUnmodeledCoefs + 1)]);
  }
}

void 
calculateProjections(SEXP regression, MERCache *cache)
{
  int *dims = DIMS_SLOT(regression);
  
  int numObservations   = dims[n_POS];
  int numUnmodeledCoefs = dims[p_POS];
  int numModeledCoefs   = dims[q_POS];
  
  double *modeledCoefProjection   = cache->modeledCoefProjection;
  double *unmodeledCoefProjection = cache->unmodeledCoefProjection;
     
  double *offDiagonalBlockRightFactorization = RZX_SLOT(regression);
  double *lowerRightBlockRightFactorization  = RX_SLOT(regression);
  
  double *denseDesignMatrix = X_SLOT(regression);
  double *response          = Y_SLOT(regression);
  
  CHM_SP rotatedSparseDesignMatrix       = A_SLOT(regression);
  CHM_FR upperLeftBlockLeftFactorization = L_SLOT(regression);
  R_CheckStack();
  
  int modelIncludesWeights = SXWT_SLOT(regression) != NULL;
  
  if (modelIncludesWeights) {
    denseDesignMatrix = cache->weightedDenseDesignMatrix;
    rotatedSparseDesignMatrix->x = Cx_SLOT(regression);
    response = cache->weightedResponse;
  } else {
    int modelOffsetsResponse = OFFSET_SLOT(regression) != NULL;
    if (modelOffsetsResponse) response = cache->weightedResponse;
  }
  
  /* solve L del1 = PAy  */ 
  // compute Perm * A * y and store in a temp
  multiplyWithPermutation(modeledCoefProjection,
                          (int *) upperLeftBlockLeftFactorization->Perm,
                          rotatedSparseDesignMatrix, response, 1);
  
  // L^-1 PAy
  solveSparseCholeskySystem(CHOLMOD_L, upperLeftBlockLeftFactorization,
                            modeledCoefProjection, 1,
                            modeledCoefProjection);
  DEBUG_PRINT_ARRAY("T~ ", modeledCoefProjection, numModeledCoefs > 10 ? 10 : numModeledCoefs);
  // solve RX' del2 = X'y - RZX'del1
  // compute X'y and store it into a temp
  multiplyMatrixByVector(denseDesignMatrix, numObservations, numUnmodeledCoefs, TRUE /* use transpose */,
                         response, unmodeledCoefProjection);
  // compute X'y - RZX' * L^-1 PAy and store back into the same temp
  multiplyMatrixByVectorWithUpdate(offDiagonalBlockRightFactorization, numModeledCoefs, numUnmodeledCoefs,
                                   TRUE /* use transpose */,
                                   modeledCoefProjection,
                                   -1.0 /* matrix product is multiplied by this before being added to the target */,
                                   unmodeledCoefProjection);
  
  if (!canProfileCommonScale(regression)) {
    Memcpy(cache->downdatedDenseResponseRotation, (const double *) unmodeledCoefProjection,
           numUnmodeledCoefs);
  }
  
  int i_one = 1;
  // result = RX^-T temp = RX^-T (X'y - RZX' L^-1 PAy)
  
  // A'x = b for A an Upper triangular, Tranposed, Non-unit matrix
  F77_CALL(dtrsv)("U", "T", "N",
                  &numUnmodeledCoefs,
                  lowerRightBlockRightFactorization,
                  &numUnmodeledCoefs,
                  unmodeledCoefProjection,
                  &i_one);
  DEBUG_PRINT_ARRAY("B~ ", unmodeledCoefProjection, numUnmodeledCoefs > 10 ? 10 : numUnmodeledCoefs);
}

void
calculatePenalizedWeightedResidualSumOfSquaresFromProjections(SEXP regression, MERCache *cache)
{
  int *dims = DIMS_SLOT(regression);
  
  int numObservations   = dims[n_POS];
  int numUnmodeledCoefs = dims[p_POS];
  int numModeledCoefs   = dims[q_POS];
  
  double *modeledCoefProjection   = cache->modeledCoefProjection;
  double *unmodeledCoefProjection = cache->unmodeledCoefProjection;
  
  int useWeightedResponse = OFFSET_SLOT(regression) || SXWT_SLOT(regression);
  double *response = (useWeightedResponse ? cache->weightedResponse : Y_SLOT(regression));
  
  cache->modeledCoefSumOfSquares = getSumOfSquares(modeledCoefProjection, numModeledCoefs);
  
  // Sadly, the penalized weighted residual sum of squares needs to be computed
  // as follows, or else, due to presumably register swaps, it does not end up
  // identical to the lmer version. And yes, there is a redundant calculation above.
  // And yes, for lmms the sum of the squares of y isn't ever going to change.
  double *deviances = DEV_SLOT(regression);
  deviances[pwrss_POS] = 
    getSumOfSquares(response, numObservations) -
    (getSumOfSquares(unmodeledCoefProjection, numUnmodeledCoefs) +
     getSumOfSquares(modeledCoefProjection, numModeledCoefs));

#ifdef PRINT_TRACE
  double ssy = getSumOfSquares(response, numObservations), sstheta = getSumOfSquares(modeledCoefProjection, numModeledCoefs), ssbeta = getSumOfSquares(unmodeledCoefProjection, numUnmodeledCoefs);
  Rprintf("ss : %llu %llu %llu %llu\n", *((unsigned long long *) (deviances + pwrss_POS)),
	      *((unsigned long long *) &ssy), *((unsigned long long *) &sstheta), *((unsigned long long *) &ssbeta));
#endif
  if (deviances[pwrss_POS] < 0.0) {
    error("Calculated PWRSS for a LMM is negative");
  }
}

void rotateProjections(SEXP regression, MERCache *cache)
{
  int *dims = DIMS_SLOT(regression);
  
  int numUnmodeledCoefs = dims[p_POS];
  int numModeledCoefs   = dims[q_POS];
  
  
  double *modeledCoef   = U_SLOT(regression);
  double *unmodeledCoef = FIXEF_SLOT(regression);
  
  Memcpy(modeledCoef,   (const double *) cache->modeledCoefProjection, numModeledCoefs);
  Memcpy(unmodeledCoef, (const double *) cache->unmodeledCoefProjection, numUnmodeledCoefs);
  
  double *lowerRightBlockRightFactorization  = RX_SLOT(regression);
  double *offDiagonalBlockRightFactorization = RZX_SLOT(regression);
  
  CHM_FR upperLeftBlockLeftFactorization = L_SLOT(regression);
  R_CheckStack();
  
  int i_one = 1;
  // solve Rx beta_tilde = beta_perp, the projection from above
  F77_CALL(dtrsv)("U", "N", "N",
                  &numUnmodeledCoefs, lowerRightBlockRightFactorization,
                  &numUnmodeledCoefs, unmodeledCoef, &i_one);
  
  // subtract from modeled coef projection RZX * beta_tilde
  multiplyMatrixByVectorWithUpdate(offDiagonalBlockRightFactorization, numModeledCoefs, numUnmodeledCoefs,
                                   FALSE /* don't use transpose */,
                                   unmodeledCoef, -1.0,
                                   modeledCoef);
  
  // u = L^-T (u - RZX beta_tilde)
  solveSparseCholeskySystem(CHOLMOD_Lt, upperLeftBlockLeftFactorization,
                            modeledCoef, 1, modeledCoef);
                            
  DEBUG_PRINT_ARRAY("T^ ", modeledCoef, numModeledCoefs > 10 ? 10 : numModeledCoefs);
  DEBUG_PRINT_ARRAY("B^ ", unmodeledCoef, numUnmodeledCoefs > 10 ? 10 : numUnmodeledCoefs);
}

void
updateRegressionForNewCommonScale(SEXP regression, MERCache *cache)
{
  int    *dims      = DIMS_SLOT(regression);
  double *deviances = DEV_SLOT(regression);
  
  int numUnmodeledCoefs = dims[p_POS];
  
  // next, we need to refactor (X'X - Rzx'Rzx + sigma^2 / sigma_beta^2 * I)
  double *lowerRightBlockRightFactorization = RX_SLOT(regression);

  // recover the cached version of X'X - LzxLzx'
  Memcpy(lowerRightBlockRightFactorization, (const double *) cache->downdatedDenseCrossproduct,
         numUnmodeledCoefs * numUnmodeledCoefs);
  
  addGaussianContributionToDenseBlock(regression, lowerRightBlockRightFactorization);
  
  getDenseCholeskyDecomposition(lowerRightBlockRightFactorization, numUnmodeledCoefs, TRIANGLE_TYPE_UPPER);

  deviances[ldRX2_POS] = 0.0;
  for (int j = 0; j < numUnmodeledCoefs; ++j) {
    deviances[ldRX2_POS] += 2.0 * log(lowerRightBlockRightFactorization[j * (numUnmodeledCoefs + 1)]);
  }
  
  // at this point, we have the correct Rx stored. now we need to
  // compute the new half projection and the new sum of squares
  double *unmodeledCoefProjection = cache->unmodeledCoefProjection;
  
  // copy in (X'Y - Rzx' theta half projection); beta half projection is Rx^-1 times that
  Memcpy(unmodeledCoefProjection, (const double *) cache->downdatedDenseResponseRotation,
         numUnmodeledCoefs);
  
  
  int i_one = 1;
  // A'x = b for A an Upper triangular, Tranposed, Non-unit matrix
  F77_CALL(dtrsv)("U", "T", "N",
                  &numUnmodeledCoefs,
                  lowerRightBlockRightFactorization,
                  &numUnmodeledCoefs,
                  unmodeledCoefProjection,
                  &i_one);

  double unmodeledCoefSumOfSquares = getSumOfSquares(unmodeledCoefProjection, numUnmodeledCoefs);
  
  // now update the penalized, weighted residual sum of squares
  deviances[pwrss_POS] = cache->responseSumOfSquares - (cache->modeledCoefSumOfSquares + unmodeledCoefSumOfSquares);
}

// As Newton's method is x_n+1 = x_n - f'(x) / f"(x), this function
// computes f'(x) and f"(x) as a function of the common scale
//
// the calculations it uses are in the accompanying pdf, but briefly
// the first derivative is related to the sample size, the residual sum of
// squares, and a new term which involves rotating the projection
// of the unmodeled coefficients
//
// the second derivative involves all of the above, plus the projection
// of the unmodeled coefficients rotated twice
//
// that is
//   f'(x) = a * N + b * PWRSS + c * || Rx^-1 beta.tilde ||^2
//   f"(x) = d * N + e * PWRSS + f * || Rx^-1 beta.tilde ||^2 + g * || Rx^-T Rx^-1 beta.tilde ||^2
static void
getDerivatives(SEXP regression, MERCache *cache,
               double *firstDerivative, double *secondDerivative)
{
  int    *dims      = DIMS_SLOT(regression);
  double *deviances = DEV_SLOT(regression);
  
  int i_one = 1;
  double d_one = 1.0;
  
  int numUnmodeledCoefs = dims[p_POS];
  
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  double *hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(unmodeledCoefPrior) + 1; // skip over the log det of the covar, not needed here
  int numHyperparameters = LENGTH(GET_SLOT(unmodeledCoefPrior, blme_prior_hyperparametersSym)) - 1;
  
  
  // take Rx and get Rx^-1
  const double *lowerRightBlockRightFactorization = RX_SLOT(regression);
  double rightFactorInverse[numUnmodeledCoefs * numUnmodeledCoefs]; // Rx^-1
  invertUpperTriangularMatrix(lowerRightBlockRightFactorization, numUnmodeledCoefs, rightFactorInverse);
  
  // calculate Lbeta^-1 * Rx^-1
  int factorIsTriangular = TRUE;
  if (numHyperparameters == 1) {
    // multiply by a scalar
    //printMatrix(lowerRightBlockRightFactorization, numUnmodeledCoefs, numUnmodeledCoefs);
    for (int col = 0; col < numUnmodeledCoefs; ++col) {
      int offset = col * numUnmodeledCoefs;
      for (int row = 0; row <= col; ++row) {
        rightFactorInverse[offset++] *= hyperparameters[0];
      }
    }
    //printMatrix(rightFactorInverse, numUnmodeledCoefs, numUnmodeledCoefs);
  } else if (numHyperparameters == numUnmodeledCoefs) {
    // left multiply by a diagonal matrix
    double *diagonal = hyperparameters;
    
    for (int col = 0; col < numUnmodeledCoefs; ++col) {
      int offset = col * numUnmodeledCoefs;
      for (int row = 0; row <= col; ++row) {
        rightFactorInverse[offset++] *= diagonal[row];
      }
    }
  } else {
    double *priorLeftFactorInverse = hyperparameters;
    // want L * R
    // Left multiply, Lower triangluar matrix, No-transpose, Non-unit
    F77_CALL(dtrmm)("L", "L", "N", "N", &numUnmodeledCoefs, &numUnmodeledCoefs, &d_one,
                    priorLeftFactorInverse, &numUnmodeledCoefs,
                    rightFactorInverse, &numUnmodeledCoefs);
    factorIsTriangular = FALSE;
  }
  
  double projectionRotation[numUnmodeledCoefs];
  Memcpy(projectionRotation, (const double *) cache->unmodeledCoefProjection, numUnmodeledCoefs);
  
  // this step corresponds to Rx^-1 * modeled coef projection
  if (factorIsTriangular) {
    // x <- Ax, A triangular. 
    // Here, A is Upper triangular, Non-tranposed, and Non-unit-diagonal

    F77_CALL(dtrmv)("U", "N", "N",
                    &numUnmodeledCoefs,
                    rightFactorInverse, &numUnmodeledCoefs,
                    projectionRotation, &i_one);
  } else {
    multiplyMatrixByVector(rightFactorInverse, numUnmodeledCoefs, numUnmodeledCoefs, FALSE,
                           projectionRotation, projectionRotation);
  }
  
  double firstRotationSumOfSquares = getSumOfSquares(projectionRotation, numUnmodeledCoefs);
  
  // now for Rx^-T Rx^-1 * modeled coef projection
  if (factorIsTriangular) {
    F77_CALL(dtrmv)("U", "T", "N",
                    &numUnmodeledCoefs,
                    rightFactorInverse, &numUnmodeledCoefs,
                    projectionRotation, &i_one);
  } else {
    multiplyMatrixByVector(rightFactorInverse, numUnmodeledCoefs, numUnmodeledCoefs, TRUE,
                           projectionRotation, projectionRotation);
  }
  
  double secondRotationSumOfSquares = getSumOfSquares(projectionRotation, numUnmodeledCoefs);
  
  
  double degreesOfFreedom = (double) (dims[isREML_POS] ? (dims[n_POS] - dims[p_POS]) : dims[n_POS]);
  double currCommonSd = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  double currCommonVariance = currCommonSd * currCommonSd;
    
  *firstDerivative =
    (deviances[pwrss_POS] / currCommonVariance - firstRotationSumOfSquares -
     degreesOfFreedom) / currCommonSd;
  *secondDerivative =
    (-3.0 * deviances[pwrss_POS] / currCommonVariance + 3.0 * firstRotationSumOfSquares +
     degreesOfFreedom) / currCommonVariance + 4.0 * secondRotationSumOfSquares;
  
  // From here, done unless REML. REML involves taking the derivative of
  // the log determinant of LxLx' (with some unmodeled covariance terms),
  // which is just the trace of the product. The second derivative is
  // the trace of the "square" of that product.
  
  if (dims[isREML_POS]) {
    int covarianceMatrixLength = numUnmodeledCoefs * numUnmodeledCoefs;
    double crossproduct[covarianceMatrixLength];
    
    // we square the left factor Lx^-T * Lx^-1. the trace of this is immediately
    // useful, but we also need the trace of its square. Fortunately, the trace
    // of AA' is simply the sum of the squares of all of the elements.
    if (factorIsTriangular) {
      singleTriangularMatrixCrossproduct(rightFactorInverse, numUnmodeledCoefs, TRUE,
                                         TRIANGLE_TYPE_UPPER, crossproduct);
    } else {
      singleMatrixCrossproduct(rightFactorInverse, numUnmodeledCoefs, numUnmodeledCoefs,
                               crossproduct, TRUE, TRIANGLE_TYPE_UPPER);
    }
    double firstOrderTrace = 0.0;
    double secondOrderTrace = 0.0;
    int offset;
    // as the cross product is symmetric, we only have to use its upper
    // triangle and the diagonal
    for (int col = 0; col < numUnmodeledCoefs; ++col) {
      offset = col * numUnmodeledCoefs;
      for (int row = 0; row < col; ++row) {
        secondOrderTrace += 2.0 * crossproduct[offset] * crossproduct[offset];
        ++offset;
      }
      
      firstOrderTrace  += crossproduct[offset];
      secondOrderTrace += crossproduct[offset] * crossproduct[offset];
    }
    
    *firstDerivative  -=  2.0 * currCommonSd * firstOrderTrace;
    *secondDerivative += -2.0 * firstOrderTrace + 4.0 * currCommonVariance * secondOrderTrace;
  }
}

// Externally callable. Sets up and fills a cache, then returns the computed
// derivatives.
//
// Works on the **deviance** scale, so -2 * logLik.
SEXP bmer_getDerivatives(SEXP regression) {
  MERCache *cache = createLMMCache(regression);

  rotateSparseDesignMatrix(regression);
  updateWeights(regression, cache);
  updateAugmentedDesignMatrixFactorizations(regression, cache);
  
  calculateProjections(regression, cache);
  calculatePenalizedWeightedResidualSumOfSquaresFromProjections(regression, cache);
   
  double firstDerivative;
  double secondDerivative;
  getDerivatives(regression, cache, &firstDerivative, &secondDerivative);
  
  deleteLMMCache(cache);
  SEXP resultExp = PROTECT(allocVector(REALSXP, 2));
  double *result = REAL(resultExp);
  result[0] = -2.0 * firstDerivative;
  result[1] = -2.0 * secondDerivative;
  UNPROTECT(1);
  
  return(resultExp);
}
