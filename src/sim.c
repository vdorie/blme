#include "sim.h"

#include <Rmath.h>
#include <Matrix.h>       // cholesky matrix types
#include <R_ext/Lapack.h>

#include "Syms.h"
#include "lmm.h"
#include "lmer.h"
#include "lmer_common.h"
#include "blmer.h"
#include "util.h"
#include "unmodeledCoefficientPrior.h"
#include "matrix.h"

#include "lmm_objectiveFunction.h"
#include "lmm_commonScale.h"

#define COMMON_SCALE_LIST_NAME          "sigma"
#define UNMODELED_COEFFICIENT_LIST_NAME "fixef"
#define MODELED_COEFFICIENT_LIST_NAME   "ranef"

static int allocateStorage(SEXP* resultPtr, SEXP regression, int numSims);
static int unmodeledCoefficientsDependOnCommonScale(SEXP regression);

static void sampleCommonScale(SEXP regression, int numSims, SEXP sims, MERCache* cache);
static void sampleUnmodeledCoefficients(SEXP regression, int numSims, SEXP sims, MERCache* cache);
static void sampleModeledCoefficients(SEXP regression, int numSims, SEXP sims);
static void addMeansToSamples(SEXP regression, int numSims, SEXP sims);

SEXP bmer_sim(SEXP regression, SEXP numSimsExp)
{
  GetRNGstate();
  
  int numSims = INTEGER(numSimsExp)[0];
  
  SEXP sims;
  int protectCount = allocateStorage(&sims, regression, numSims);
  
  MERCache* cache = NULL;
  if (unmodeledCoefficientsDependOnCommonScale(regression)) {
    cache = createLMMCache(regression);
    calculateFirstHalfOfProjections(regression, cache);
  }
  
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  
  if (isLinearModel) sampleCommonScale(regression, numSims, sims, cache);
  
  // samples them as mean 0
  sampleUnmodeledCoefficients(regression, numSims, sims, cache);
  sampleModeledCoefficients  (regression, numSims, sims);
  
  addMeansToSamples(regression, numSims, sims);
  
  if (cache != NULL) deleteLMMCache(cache);
  
  UNPROTECT(protectCount);
  
  PutRNGstate();
  
  return(sims);
}

static void sampleCommonScale(SEXP regression, int numSims, SEXP sims, MERCache* cache)
{
  int    *dims      = DIMS_SLOT(regression);
  double *deviances = DEV_SLOT(regression);
  
  double* commonScaleSamples = REAL(getListElement(sims, COMMON_SCALE_LIST_NAME));
  
  enum { SIM_TYPE_DIRECT = 0, SIM_TYPE_FIXED, SIM_TYPE_APPROX  } simType = SIM_TYPE_DIRECT;
  
  double degreesOfFreedom = (double) (dims[n_POS] - dims[p_POS]);
  
  SEXP commonScalePrior = GET_SLOT(regression, blme_commonScalePriorSym);
  priorType_t commonScalePriorType = PRIOR_TYPE_SLOT(commonScalePrior);
  
  if (commonScalePriorType == PRIOR_TYPE_DIRECT) {
    priorFamily_t priorFamily = PRIOR_FAMILIES_SLOT(commonScalePrior)[0];
    
    if (priorFamily == PRIOR_FAMILY_POINT) simType = SIM_TYPE_FIXED;
    else if (priorFamily == PRIOR_FAMILY_INVGAMMA) {
      simType = SIM_TYPE_DIRECT; // possibly, is conjugate, but unmodeled coef prior can muck that up
      degreesOfFreedom += 2.0 * (PRIOR_HYPERPARAMETERS_SLOT(commonScalePrior)[0] + 1.0);
    }
  }
  
  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  priorType_t unmodeledCoefficientPriorType = PRIOR_TYPE_SLOT(unmodeledCoefficientPrior);
  
  if (simType != SIM_TYPE_FIXED &&
      unmodeledCoefficientPriorType == PRIOR_TYPE_DIRECT &&
      PRIOR_FAMILIES_SLOT(unmodeledCoefficientPrior)[0] == PRIOR_FAMILY_GAUSSIAN)
  {
    int priorScale = PRIOR_SCALES_SLOT(unmodeledCoefficientPrior)[0];
    priorCommonScale_t onCommonScale = getCommonScaleBit(priorScale);
    if (onCommonScale) {
      // prior on common scale gets extra dof
      degreesOfFreedom += (double) dims[p_POS]; 
    } else {
      // prior on absolute scale means that marginal posterior requires approximation
      simType = SIM_TYPE_APPROX;
    }
  }
  
  switch (simType) {
    case SIM_TYPE_FIXED:
    {
      double commonScale = (dims[isREML_POS] ? deviances[sigmaREML_POS] : deviances[sigmaML_POS]);
      for (int i = 0; i < numSims; ++i) commonScaleSamples[i] = commonScale;
      break;
    } 
    case SIM_TYPE_APPROX:
    {
      double commonScale = (dims[isREML_POS] ? deviances[sigmaREML_POS] : deviances[sigmaML_POS]);
      double commonVariance = commonScale * commonScale;
      
      // do inv-gamma matching posterior in mode and curvature at mode
      double firstDerivative, secondDerivative;
      getCommonScaleDerivatives(regression, cache, &firstDerivative, &secondDerivative);
      
      // the mess is because the derivatives are of sigma, not sigma.sq, sigma.sq being
      // the scale on which the posterior looks like an inv-gamma.
      //
      // Admittedly, the first derivative should be 0. Numerically, it's close, but 'eh.
      double curvature = (secondDerivative - firstDerivative / commonScale) / (4.0 * commonVariance);
      
      double shape = -commonVariance * commonVariance * curvature;
      double scale = 1.0 / (shape * commonVariance);
      shape -= 1.0;
      
      // who would have thought that an undocumented function treats its parameters in a fashion
      // opposite that of "default" R?. apparently, rgamma samples propto x^(a - 1) exp(-x/b), despite,
      // in R: rgamma(1, a, b) propto x^(a - 1) exp(-x*b). Thanks R!
      for (int i = 0; i < numSims; ++i) commonScaleSamples[i] = sqrt(1.0 / rgamma(shape, scale));
      
      break;
    }
    default:
    {
      double shape = degreesOfFreedom / 2.0;
      double scale = 1.0 / (deviances[pwrss_POS] / 2.0);
      for (int i = 0; i < numSims; ++i) commonScaleSamples[i] = sqrt(1.0 / rgamma(shape, scale));
      break;
    }
  }
}

static void sampleUnmodeledCoefficients(SEXP regression, int numSims, SEXP sims, MERCache* cache)
{
  int    *dims      = DIMS_SLOT(regression);
  
  int numUnmodeledCoefs = dims[p_POS];
  
  double* unmodeledCoefficientSamples = REAL(getListElement(sims, UNMODELED_COEFFICIENT_LIST_NAME));
  double* commonScaleSamples          = REAL(getListElement(sims, COMMON_SCALE_LIST_NAME));
  
  if (unmodeledCoefficientsDependOnCommonScale(regression)) {
    // we cache X'X - Rzx'Rzx
    // chol factor that we want is X'X - Rzx'Rzx + sigma^2 Sigma.beta^-1
    int blockMatrixSize = numUnmodeledCoefs * numUnmodeledCoefs;
    double* lowerRightBlock     = Alloca(blockMatrixSize, double);
    double* baseLowerRightBlock = Alloca(blockMatrixSize, double); 
    
    computeDowndatedDenseCrossproduct(regression, cache, baseLowerRightBlock);
    
    long long offset = 0;
    int i_one = 1;
    for (int i = 0; i < numSims; ++i) {
      for (int j = 0; j < numUnmodeledCoefs; ++j) {
        unmodeledCoefficientSamples[offset++] = norm_rand() * commonScaleSamples[i];
      }
      
      
      Memcpy(lowerRightBlock, (double* const) baseLowerRightBlock, blockMatrixSize);
      
      // prior contribution
      addGaussianContributionToDenseBlock(regression, lowerRightBlock, commonScaleSamples[i]);
      
      // factorized
      int choleskyResult = getDenseCholeskyDecomposition(lowerRightBlock, numUnmodeledCoefs, TRIANGLE_TYPE_UPPER);
      if (choleskyResult > 0) error("Leading minor %d of downdated X'X is not positive definite.", choleskyResult);
      if (choleskyResult < 0) error("Illegal argument %d to cholesky decomposition (dpotrf).", -choleskyResult);
      
      // multiply sims by left factor
      F77_CALL(dtrsv)("U", "T", "N", &numUnmodeledCoefs, lowerRightBlock,
                      &numUnmodeledCoefs, unmodeledCoefficientSamples + i * numUnmodeledCoefs, &i_one);
    }
  } else {
    double* lowerRightBlock = RX_SLOT(regression);
    
    long long offset = 0;
    for (int i = 0; i < numSims; ++i) {
      for (int j = 0; j < numUnmodeledCoefs; ++j) {
        unmodeledCoefficientSamples[offset++] = norm_rand() * commonScaleSamples[i];
      }
    }
    
    double d_one = 1.0;
    // multiply sims by left factor
    F77_CALL(dtrsm)("L", "U", "N", "N", &numUnmodeledCoefs, &numSims, &d_one, lowerRightBlock, &numUnmodeledCoefs,
                    unmodeledCoefficientSamples, &numUnmodeledCoefs);
  }
}

static void sampleModeledCoefficients(SEXP regression, int numSims, SEXP sims)
{
  int* dims = DIMS_SLOT(regression);
  int numModeledCoefs   = dims[q_POS];
  int numUnmodeledCoefs = dims[p_POS];
  
  double* commonScaleSamples          = REAL(getListElement(sims, COMMON_SCALE_LIST_NAME));
  double* unmodeledCoefficientSamples = REAL(getListElement(sims, UNMODELED_COEFFICIENT_LIST_NAME));
  double*   modeledCoefficientSamples = REAL(getListElement(sims,   MODELED_COEFFICIENT_LIST_NAME));
  
  double* upperRightFactor = RZX_SLOT(regression);
  
  long long offset = 0;
  for (int i = 0; i < numSims; ++i) {
    for (int j = 0; j < numModeledCoefs; ++j) {
      modeledCoefficientSamples[offset++] = norm_rand() * commonScaleSamples[i];
    }
  }
  
  // if we have samples of fixef, ranef samples are noise + a transformation of those samples
  //
  // ranef samp = Lz^-T (noise - Rzx fixef samp)
  
  // forms noise - Rzx fixef sampl
  double d_one = 1.0, d_minusOne = -1.0;
  // C:= alpha op(A) * op(B) + beta * C
  F77_CALL(dgemm)("No transpose, A op", "No transpose, B op",
                  &numModeledCoefs /* M = num rows A */ , &numSims /* N = num cols B */,
                  &numUnmodeledCoefs /* K = num cols A */, &d_one /* alpha */,
                  upperRightFactor /* A */, &numModeledCoefs /* "stride" of A */ ,
                  unmodeledCoefficientSamples /* B */, &numUnmodeledCoefs /* stride of B */,
                  &d_minusOne /* beta */,
                  modeledCoefficientSamples /* C */, &numModeledCoefs /* stride of C */);
  
  
  CHM_FR upperLeftFactor = L_SLOT(regression); // allocates on stack
  R_CheckStack();
  // now invert and multiply by Lz'
  solveSparseCholeskySystem(CHOLMOD_Lt, upperLeftFactor,
                            modeledCoefficientSamples /* right hand side */, numSims,
                            modeledCoefficientSamples /* target */);
  
  // at this point, we have spherical ranef simulations
  // we need to scale by their covariance (and reverse the fill-reducing Cholmod permutation)
  
  
  int numFactors = dims[nt_POS];
  
  SparseMatrixStructure sparseStructure;
  sparseStructure.factorDimensions = Alloca(numFactors, int);
  sparseStructure.numGroupsPerFactor = Alloca(numFactors, int);
  double** stMatrices = Alloca(numFactors, double*);
  double* tempColumn = Alloca(numModeledCoefs, double);
  R_CheckStack();
  
  int* sparseRowForFactor = Gp_SLOT(regression);
  getSparseContentAndStructure(GET_SLOT(regression, lme4_STSym), sparseRowForFactor,
                               stMatrices, &sparseStructure);
  int* cholmodPerm = PERM_VEC(regression);
  
  double* modeledCoefficientColumn = modeledCoefficientSamples;
  for (int i = 0; i < numSims; ++i) {
    
    // permutation first
    Memcpy(tempColumn, (double* const) modeledCoefficientColumn, numModeledCoefs);
    for (int j = 0; j < numModeledCoefs; ++j) {
      modeledCoefficientColumn[cholmodPerm[j]] = tempColumn[j];
    }
    
    
    // now to rotate/rescale by left factor = T*S
    for (int j = 0; j < numFactors; ++j) {
      int factorDimension = sparseStructure.factorDimensions[j];
      int numGroups       = sparseStructure.numGroupsPerFactor[j];
      
      for (int k = 0; k < factorDimension; ++k) {
        
        // multiply by S
        double scale = stMatrices[j][k * (factorDimension + 1)]; // pull scale off of diagonal
        int base = sparseRowForFactor[j] + k * numGroups;
        for (int l = 0; l < numGroups; ++l) {
          modeledCoefficientColumn[base + l] *= scale;
        }
        
        // multiply by T; this is slightly hacky, in that it takes the column vector
        // and pretends its a matrix. Each "column" corresponds to the replications of a single
        // ranef, and each row a single group in that factor e.g.:
        // 
        // for factor 1:
        //   group1:  int slo1 slo2
        //   group2:  int slo1 slo2
        //   group3:  int slo1 slo2
        //
        // the way T is stored corresponds to:
        //   g1:      int
        //            slo1
        //            slo2
        //   g2:      int
        //            slo1
        //            slo2
        //
        // so that, yeah, the ranef vector, matrixified, is just the transpose of the form that
        // T is in
        if (factorDimension > 1) {
          // B:= alpha * B * op(A)
          F77_CALL(dtrmm)("Right side mult", "Lower triangular", "Transpose", "Unit triangular",
                          &numGroups /* M = numRows B */, &factorDimension /* N = numCols B */, &d_one /* alpha */,
                          stMatrices[j] /* A */, &factorDimension /* stride of A */,
                          modeledCoefficientColumn + sparseRowForFactor[j] /* B */, &numGroups /* stride of B */);
          
        }
      }
    }
        
    modeledCoefficientColumn += numModeledCoefs;
  }
}

// adds fixef and ranef to the mean 0 sims we have already generated
static void addMeansToSamples(SEXP regression, int numSims, SEXP sims)
{
  int* dims = DIMS_SLOT(regression);
  int numModeledCoefs   = dims[q_POS];
  int numUnmodeledCoefs = dims[p_POS];
  
  double* unmodeledCoefficientSamples = REAL(getListElement(sims, UNMODELED_COEFFICIENT_LIST_NAME));
  double*   modeledCoefficientSamples = REAL(getListElement(sims,   MODELED_COEFFICIENT_LIST_NAME));
  
  double* unmodeledCoefficients = FIXEF_SLOT(regression);
  double*   modeledCoefficients = RANEF_SLOT(regression);
  
  double* unmodeledCoefficientColumn = unmodeledCoefficientSamples;
  double*   modeledCoefficientColumn =   modeledCoefficientSamples;
  for (int i = 0; i < numSims; ++i) {
    for (int j = 0; j < numUnmodeledCoefs; ++j) {
      unmodeledCoefficientColumn[j] += unmodeledCoefficients[j];
    }
    for (int j = 0; j < numModeledCoefs; ++j) {
      modeledCoefficientColumn[j] += modeledCoefficients[j];
    }
    
    unmodeledCoefficientColumn += numUnmodeledCoefs;
    modeledCoefficientColumn   += numModeledCoefs;
  }
}

static int
unmodeledCoefficientsDependOnCommonScale(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  if (!isLinearModel) return 0;
  
  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  priorType_t unmodeledCoefficientPriorType = PRIOR_TYPE_SLOT(unmodeledCoefficientPrior);
  
  return (unmodeledCoefficientPriorType == PRIOR_TYPE_DIRECT &&
          PRIOR_FAMILIES_SLOT(unmodeledCoefficientPrior)[0] == PRIOR_FAMILY_GAUSSIAN &&
          getCommonScaleBit(PRIOR_SCALES_SLOT(unmodeledCoefficientPrior)[0]) == PRIOR_COMMON_SCALE_FALSE);
}

static
int allocateStorage(SEXP* resultPtr, SEXP regression, int numSims)
{
  int protectCount = 0;
  PROTECT(*resultPtr = allocVector(VECSXP, 3));
  ++protectCount;
  
  SEXP result = *resultPtr;
  
  int* dims = DIMS_SLOT(regression);
  
  SEXP fixef, ranef, sigma;
  PROTECT(fixef = allocVector(REALSXP, numSims * dims[p_POS]));
  ++protectCount;
  SET_DIMS(fixef, dims[p_POS], numSims);
  
  
  PROTECT(ranef = allocVector(REALSXP, numSims * dims[q_POS]));
  ++protectCount;
  SET_DIMS(ranef, dims[q_POS], numSims);
  
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  if (isLinearModel) {
    PROTECT(sigma = allocVector(REALSXP, numSims));
    ++protectCount;
  } else {
    sigma = ScalarReal(NA_REAL);
  }
  
  SEXP resultNames;
  setAttrib(result, R_NamesSymbol, resultNames = allocVector(STRSXP, 3));
  
  SET_VECTOR_ELT(result, 0, fixef);
  SET_STRING_ELT(resultNames, 0, mkChar(UNMODELED_COEFFICIENT_LIST_NAME));
  SET_VECTOR_ELT(result, 1, ranef);
  SET_STRING_ELT(resultNames, 1, mkChar(MODELED_COEFFICIENT_LIST_NAME));
  SET_VECTOR_ELT(result, 2, sigma);
  SET_STRING_ELT(resultNames, 2, mkChar(COMMON_SCALE_LIST_NAME));
  
  return protectCount;
}
