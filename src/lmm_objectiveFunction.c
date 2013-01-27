#include "lmm.h"
#include "lmm_objectiveFunction.h"

#include <Rmath.h>
#include <Matrix.h>
#include <R_ext/Lapack.h>

#include "__lmmMerCache.h"

#include "Syms.h"
#include "lmer.h"
#include "blmer.h"

#include "lmer_common.h"
#include "util.h"
#include "matrix.h"
#include "common_inlines.h"

#include "unmodeledCoefficientPrior.h"
#include "lmm_commonScale.h"
#include "lmm_priors.h"

extern	       // cholmod_common struct initialized in R_init_lme4 
cholmod_common cholmodCommon;

// not static so that they can be called by test code
void calculateObjectiveFunction(SEXP regression, MERCache* cache, const double* parameters);
void calculateDeviances(SEXP regression, MERCache* cache);

double lmmGetObjectiveFunction(SEXP regression, MERCache* cache, const double* parameters)
{
  weightMatrices(regression, cache);
  updateMatrixFactorizations(regression, cache);
  setPriorVaryingContributionsToCommonScale(regression, cache, parameters);
  
  // these morally speaking calculate the mode of the joint distribution
  calculateFirstHalfOfProjections(regression, cache); // contains the values needed to profile out the coefficients
  optimizeCommonScale(regression, cache);
    
  calculateObjectiveFunction(regression, cache, parameters);
  
  return(cache->objectiveFunctionValue);
}

double lmmGetObjectiveFunctionForFixedCommonScale(SEXP regression, MERCache* cache, const double* parameters)
{
  weightMatrices(regression, cache);
  updateMatrixFactorizations(regression, cache);
  setPriorVaryingContributionsToCommonScale(regression, cache, parameters);
  
  calculateFirstHalfOfProjections(regression, cache);
  
  calculateObjectiveFunction(regression, cache, parameters);
  
  return(cache->objectiveFunctionValue);
}


void lmmCalculateProjectionsForSingleArgumentAnova(SEXP regression,
                                                   double *modeledCoefProjection,
                                                   double *unmodeledCoefProjection)
{
  MERCache* cache = createLMMCache(regression);
  
  // replace where the cache points to save a little bit of time
  // as opposed to copying from the temp
  double* cacheModeledCoefProjection   = cache->modeledCoefProjection;
  double* cacheUnmodeledCoefProjection = cache->unmodeledCoefProjection;
  
  cache->modeledCoefProjection   = modeledCoefProjection;
  cache->unmodeledCoefProjection = unmodeledCoefProjection;
  
  weightMatrices(regression, cache);
  updateMatrixFactorizations(regression, cache);
  calculateFirstHalfOfProjections(regression, cache);
  
  // since they'll be freed when we delete
  cache->modeledCoefProjection   = cacheModeledCoefProjection;
  cache->unmodeledCoefProjection = cacheUnmodeledCoefProjection;
  
  deleteLMMCache(cache);
}

void lmmFinalizeOptimization(SEXP regression, MERCache* cache)
{
  calculateSecondHalfOfProjections(regression, cache);
  calculateSumsOfSquares(regression, cache);
  calculateDeviances(regression, cache);
}

double lmmGetPriorPenalty(SEXP regression, MERCache* cache, const double* parameters)
{
  if (!ISNAN(cache->priorContribution)) return(cache->priorContribution);

  // fall back to blmer function default, since cache hasn't been setup
  return(cache->priorContribution = getPriorPenalty(regression, cache, parameters));
}

// forward decls used by weightMatrices
//static void updateResidualWeights(SEXP regression);
//static void updateObservationWeights(SEXP regression);
static void weightSparseDesignMatrix(SEXP regression);
static void weightResponseAndDenseDesignMatrix(SEXP regression, MERCache *cache);

void weightMatrices(SEXP regression, MERCache* cache)
{
  // see function definitions as to why these are commented out
  // updateResidualWeights(regression);
  // updateObservationWeights(regression);
  
  weightSparseDesignMatrix(regression);
  weightResponseAndDenseDesignMatrix(regression, cache);
}


// used by updateMatrixFactorizations
static void updateUpperLeftCholeskyBlock(SEXP regression);
static void updateOffDiagonalCholeskyBlock(SEXP regression, MERCache *cache);
static void updateLowerRightCholeskyBlock(SEXP regression, MERCache *cache);

void updateMatrixFactorizations(SEXP regression, MERCache* cache)
{
  updateUpperLeftCholeskyBlock(regression);
  
#ifdef PRINT_TRACE
  const int* dims = DIMS_SLOT(regression);
  CHM_FR L = L_SLOT(regression); R_CheckStack();
#endif
  
  DEBUG_PRINT_ARRAY("L.x", L->x, ((int *) L->p)[dims[q_POS]] > 10 ? 10 : ((int *) L->p)[dims[q_POS]]);
  
  updateOffDiagonalCholeskyBlock(regression, cache);
  
  DEBUG_PRINT_ARRAY("RZX", RZX_SLOT(regression), dims[p_POS] * dims[q_POS] > 10 ? 10 : dims[p_POS] * dims[q_POS]);
  
  updateLowerRightCholeskyBlock(regression, cache);
  
  DEBUG_PRINT_ARRAY("RX ", RX_SLOT(regression), dims[p_POS] * dims[p_POS] > 10 ? 10 : dims[p_POS] * dims[p_POS]);
}


void calculateFirstHalfOfProjections(SEXP regression, MERCache* cache)
{
  const int* dims = DIMS_SLOT(regression);
  
  int numObservations   = dims[n_POS];
  int numUnmodeledCoefs = dims[p_POS];
  int numModeledCoefs   = dims[q_POS];
  
  double* modeledCoefProjection   = cache->modeledCoefProjection;
  double* unmodeledCoefProjection = cache->unmodeledCoefProjection;
  
  double* offDiagonalBlockRightFactorization = RZX_SLOT(regression);
  double* lowerRightBlockRightFactorization  = RX_SLOT(regression);
  
  double* denseDesignMatrix = X_SLOT(regression);
  double* response          = Y_SLOT(regression);
  
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
  applyMatrixToVector(denseDesignMatrix, numObservations, numUnmodeledCoefs, TRUE /* use transpose */,
                      response, unmodeledCoefProjection);
  // compute X'y - RZX' * L^-1 PAy and store back into the same temp
  applyMatrixToVectorWithUpdate(offDiagonalBlockRightFactorization, numModeledCoefs, numUnmodeledCoefs,
                                TRUE /* use transpose */,
                                modeledCoefProjection,
                                -1.0 /* matrix product is multiplied by this before being added to the target */,
                                unmodeledCoefProjection);
  
  if (cache->commonScaleOptimization == CSOT_BRUTE_FORCE) {
    // with a brute force optimization, we'll potentially be re-using this bit and changing the
    // next by multiplying in the common scale
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
  
  cache->modeledCoefProjectionSumOfSquares   = getSumOfSquares(modeledCoefProjection, numModeledCoefs);
  cache->unmodeledCoefProjectionSumOfSquares = getSumOfSquares(unmodeledCoefProjection, numUnmodeledCoefs);
  cache->totalSumOfSquares = cache->responseSumOfSquares -
     (cache->modeledCoefProjectionSumOfSquares + cache->unmodeledCoefProjectionSumOfSquares);
}

void calculateObjectiveFunction(SEXP regression, MERCache* cache, const double* parameters)
{
  const int*    dims      = DIMS_SLOT(regression);
  const double* deviances = DEV_SLOT(regression);
  
  double degreesOfFreedom = (double) (dims[n_POS] - (dims[isREML_POS] ? dims[p_POS] : 0.0));
  
  // constants
  double result = cache->priorDevianceConstantPart + 2.0 * degreesOfFreedom * M_LN_SQRT_2PI;
  
  // determinants
  result += deviances[ldL2_POS] + (dims[isREML_POS] ? deviances[ldRX2_POS] : 0.0);
  
  double sigma = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  double sigma_sq = sigma * sigma;
  
  degreesOfFreedom += cache->priorCommonScaleDegreesOfFreedom;
  
  // common scale polynomial term
  result += degreesOfFreedom * log(sigma_sq);
  
  // exponential term related to sigma
  double exponentialTerm = (cache->totalSumOfSquares +
                            2.0 * (cache->mTwoExponentialTermConstantPart + cache->mTwoExponentialTermVaryingPart)) / sigma_sq;
  exponentialTerm += 2.0 * (cache->mOneExponentialTermConstantPart + cache->mOneExponentialTermVaryingPart) / sigma;
  exponentialTerm += 2.0 * (cache->oneExponentialTermConstantPart + cache->oneExponentialTermVaryingPart) * sigma;
  exponentialTerm += 2.0 * (cache->twoExponentialTermConstantPart + cache->twoExponentialTermVaryingPart) * sigma_sq;
  
  result += exponentialTerm;
  
  // any lingering covariance prior terms
  result += getPriorDevianceCommonScaleFreeVaryingPart(regression, parameters); 
  
  cache->objectiveFunctionValue = result;
}

void calculateSecondHalfOfProjections(SEXP regression, MERCache *cache)
{
  const int* dims = DIMS_SLOT(regression);
  
  int numUnmodeledCoefs = dims[p_POS];
  int numModeledCoefs   = dims[q_POS];
  
  
  double* modeledCoef   = U_SLOT(regression);
  double* unmodeledCoef = FIXEF_SLOT(regression);
  
  Memcpy(modeledCoef,   (const double*) cache->modeledCoefProjection, numModeledCoefs);
  Memcpy(unmodeledCoef, (const double*) cache->unmodeledCoefProjection, numUnmodeledCoefs);
  
  double* lowerRightBlockRightFactorization  = RX_SLOT(regression);
  double* offDiagonalBlockRightFactorization = RZX_SLOT(regression);
  
  CHM_FR upperLeftBlockLeftFactorization = L_SLOT(regression);
  R_CheckStack();
  
  int i_one = 1;
  // solve Rx beta_tilde = beta_perp, the projection from above
  F77_CALL(dtrsv)("U", "N", "N",
                  &numUnmodeledCoefs, lowerRightBlockRightFactorization,
                  &numUnmodeledCoefs, unmodeledCoef, &i_one);
  
  // subtract from modeled coef projection RZX * beta_tilde
  applyMatrixToVectorWithUpdate(offDiagonalBlockRightFactorization, numModeledCoefs, numUnmodeledCoefs,
                                FALSE /* don't use transpose */,
                                unmodeledCoef, -1.0,
                                modeledCoef);
  
  // u = L^-T (u - RZX beta_tilde)
  solveSparseCholeskySystem(CHOLMOD_Lt, upperLeftBlockLeftFactorization,
                            modeledCoef, 1, modeledCoef);
  
  DEBUG_PRINT_ARRAY("T^ ", modeledCoef, numModeledCoefs > 10 ? 10 : numModeledCoefs);
  DEBUG_PRINT_ARRAY("B^ ", unmodeledCoef, numUnmodeledCoefs > 10 ? 10 : numUnmodeledCoefs);
}

// the notion of the relevant sums of squares differs depending on the priors
// see the accompanying linear optimization pdf
void calculateSumsOfSquares(SEXP regression, MERCache* cache)
{
  const int* dims   = DIMS_SLOT(regression);
  double* deviances = DEV_SLOT(regression);
  
  int numUnmodeledCoefs = dims[p_POS];
  int numModeledCoefs   = dims[q_POS];
  
  const double* unmodeledCoefs = FIXEF_SLOT(regression);
  
  deviances[usqr_POS] = getSumOfSquares(U_SLOT(regression), numModeledCoefs);
  
  
  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  if (PRIOR_TYPE_SLOT(unmodeledCoefficientPrior)        != PRIOR_TYPE_DIRECT ||
      PRIOR_FAMILIES_SLOT(unmodeledCoefficientPrior)[0] != PRIOR_FAMILY_GAUSSIAN)
  {
    deviances[pwrss_POS] = cache->totalSumOfSquares;
    deviances[wrss_POS]  = deviances[disc_POS] = deviances[pwrss_POS] - deviances[usqr_POS];
  } else {
    double unmodeledCoefficientSumOfSquares = getUnmodeledCoefficientDensityExponentialPart(unmodeledCoefficientPrior, unmodeledCoefs, numUnmodeledCoefs);
    
    if (getCommonScaleBit(PRIOR_SCALES_SLOT(unmodeledCoefficientPrior)[0]) == PRIOR_COMMON_SCALE_FALSE) {
      double commonScale = DEV_SLOT(regression)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
      commonScale *= commonScale;
      
      deviances[pwrss_POS] = cache->totalSumOfSquares - commonScale * unmodeledCoefficientSumOfSquares;
      deviances[wrss_POS]  = deviances[disc_POS] = deviances[pwrss_POS] - deviances[usqr_POS];
    } else {
      deviances[pwrss_POS] = cache->totalSumOfSquares;
      deviances[wrss_POS]  = deviances[disc_POS] = deviances[pwrss_POS] - (unmodeledCoefficientSumOfSquares + deviances[usqr_POS]);
    }
  }
  
  deviances[pwrss_POS] += 2.0 * (cache->mTwoExponentialTermConstantPart + cache->mTwoExponentialTermVaryingPart);
}

void calculateDeviances(SEXP regression, MERCache* cache)
{
  const int* dims   = DIMS_SLOT(regression);
  double* deviances = DEV_SLOT(regression);
  
  double exponentialTerm = deviances[usqr_POS] + deviances[wrss_POS];
  
  double degreesOfFreedom = (double) dims[n_POS];
  double sigma_sq = deviances[sigmaML_POS] * deviances[sigmaML_POS];
  
  deviances[  ML_POS] = deviances[ldL2_POS] + 
    degreesOfFreedom * (2.0 * M_LN_SQRT_2PI + log(sigma_sq)) + exponentialTerm / sigma_sq;
  
  degreesOfFreedom -= (double) dims[p_POS];
  sigma_sq = deviances[sigmaREML_POS] * deviances[sigmaREML_POS];
  
 deviances[REML_POS] = deviances[ldL2_POS] + deviances[ldRX2_POS] +
    degreesOfFreedom * (2.0 * M_LN_SQRT_2PI + log(sigma_sq)) + exponentialTerm / sigma_sq;
  
  cache->priorContribution = cache->objectiveFunctionValue - deviances[dims[isREML_POS] ? REML_POS : ML_POS];
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
  
  double* denseDesignMatrix         = X_SLOT(regression);
  double* weightedDenseDesignMatrix = cache->weightedDenseDesignMatrix;
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

// matrix factorization helpers
static void updateUpperLeftCholeskyBlock(SEXP regression)
{
  CHM_SP rotatedSparseDesignMatrix       = A_SLOT(regression);
  CHM_FR upperLeftBlockLeftFactorization = L_SLOT(regression);
  R_CheckStack();
  
  int modelIncludesWeights = SXWT_SLOT(regression) != NULL;
  
  // replace A with C if needed, as previously computed
  if (modelIncludesWeights) {
    double* cMatrixExplicitStorage = Cx_SLOT(regression);
    
    if (cMatrixExplicitStorage == NULL) {
      error("For an lmm with weights, Cx is required to exist but has not been defined.");
    }
    rotatedSparseDesignMatrix->x = (void*) cMatrixExplicitStorage;
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
  
  double* deviances = DEV_SLOT(regression);
  deviances[ldL2_POS] = M_chm_factor_ldetL2(upperLeftBlockLeftFactorization);
}

static void updateOffDiagonalCholeskyBlock(SEXP regression, MERCache* cache)
{
  const int* dims = DIMS_SLOT(regression);
  int numUnmodeledCoefs = dims[p_POS];
  
  double* denseDesignMatrix = X_SLOT(regression);
  double* offDiagonalBlockRightFactorization = RZX_SLOT(regression);
  
  CHM_SP rotatedSparseDesignMatrix       = A_SLOT(regression);
  CHM_FR upperLeftBlockLeftFactorization = L_SLOT(regression);
  R_CheckStack();
  
  int modelIncludesWeights = SXWT_SLOT(regression) != NULL;
  
  // replace model matrices by their weighted versions if necessary
  if (modelIncludesWeights) {
    denseDesignMatrix = cache->weightedDenseDesignMatrix;
    
    double* cMatrixExplicitStorage = Cx_SLOT(regression);
    if (cMatrixExplicitStorage == NULL) {
      error("For an lmm with weights, Cx is required to exist but has not been defined.");
    }
    
    rotatedSparseDesignMatrix->x = (void*) cMatrixExplicitStorage;
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
  const int* dims = DIMS_SLOT(regression);
  double commonScale = DEV_SLOT(regression)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  
  int numUnmodeledCoefs = dims[p_POS];
  
  double* lowerRightBlockRightFactorization  = RX_SLOT(regression);
  
  computeDowndatedDenseCrossproduct(regression, cache, lowerRightBlockRightFactorization);
  
  if (cache->commonScaleOptimization == CSOT_BRUTE_FORCE) {
    Memcpy(cache->downdatedDenseCrossproduct, (const double*) lowerRightBlockRightFactorization,
           numUnmodeledCoefs * numUnmodeledCoefs);
  }
  
  addGaussianContributionToDenseBlock(regression, lowerRightBlockRightFactorization, commonScale);
  
  // factor result
  int choleskyResult = getDenseCholeskyDecomposition(lowerRightBlockRightFactorization, numUnmodeledCoefs,
                                                     TRIANGLE_TYPE_UPPER);
  
  if (choleskyResult > 0) error("Leading minor %d of downdated X'X is not positive definite.", choleskyResult);
  if (choleskyResult < 0) error("Illegal argument %d to cholesky decomposition (dpotrf).", -choleskyResult);
  
  double* deviances = DEV_SLOT(regression);
  // accumulate log(det(RX)^2)
  // matrix is triangular
  deviances[ldRX2_POS] = 0.0;
  for (int j = 0; j < numUnmodeledCoefs; ++j) {
    deviances[ldRX2_POS] += 2.0 * log(lowerRightBlockRightFactorization[j * (numUnmodeledCoefs + 1)]);
  }
}

// populates the cache so that the derivative code can execute correctly;
// mimics a call to lmmGetObjectiveFunction up until common scale optimization
// takes place
void _getCommonScaleDerivatives(SEXP regression, MERCache* cache, const double* parameters, double* firstDerivative, double* secondDerivative)
{
  weightMatrices(regression, cache);
  updateMatrixFactorizations(regression, cache);
  setPriorVaryingContributionsToCommonScale(regression, cache, parameters);
  
  calculateFirstHalfOfProjections(regression, cache);
  
  getCommonScaleDerivatives(regression, cache, firstDerivative, secondDerivative);
}

// populates the cache so that the derivative code can execute correctly;
// mimics a call to lmmGetObjectiveFunction but doesn't store estimated
// common scale, returning it instead
double _getOptimalCommonScale(SEXP regression, MERCache* cache, const double* parameters)
{
  double* deviances = DEV_SLOT(regression);
  
  weightMatrices(regression, cache);
  updateMatrixFactorizations(regression, cache);
  setPriorVaryingContributionsToCommonScale(regression, cache, parameters);
  
  calculateFirstHalfOfProjections(regression, cache);
  
  double sigmaML_old   = deviances[sigmaML_POS];
  double sigmaREML_old = deviances[sigmaREML_POS];
  
  optimizeCommonScale(regression, cache);

  double result = deviances[DIMS_SLOT(regression)[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  deviances[sigmaML_POS]   = sigmaML_old;
  deviances[sigmaREML_POS] = sigmaREML_old;

  return(result);
}
