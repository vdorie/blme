#include "lmm_commonScale.h"

#include <R_ext/Lapack.h>

#include "__lmmMerCache.h"
#include "lmer.h"
#include "Syms.h"
#include "blmer_types.h"

#include "common_inlines.h"
#include "matrix.h"
#include "unmodeledCoefficientPrior.h"

static void getDerivativesOfSumOfSquares(SEXP regression, MERCache* cache, double* firstDerivative, double* secondDerivative);

void optimizeCommonScale(SEXP regression, MERCache* cache)
{
  const double* deviances = DEV_SLOT(regression);
  const int*    dims      = DIMS_SLOT(regression);
  
  int commonScaleRequiresBruteForce = cache->commonScaleOptimization == CSOT_BRUTE_FORCE;
  int commonScaleCanBeProfiled      = cache->commonScaleOptimization != CSOT_BRUTE_FORCE &&
                                      cache->commonScaleOptimization != CSOT_NA;
  
  // Rprintf("csot: %d, brute: %d, prof: %d\n", cache->commonScaleOptimization, commonScaleRequiresBruteForce, commonScaleCanBeProfiled);
  
  if (commonScaleRequiresBruteForce) {
    double currCommonScale = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
    double prevCommonScale;
    
    do {
      prevCommonScale = currCommonScale;
      
      // implictly updates the dense factorization and the half-projection that depends on it
      // also updates the deviances array
      currCommonScale = performOneStepOfNewtonsMethodForCommonScale(regression, cache);
    } while (fabs(currCommonScale - prevCommonScale) >= COMMON_SCALE_OPTIMIZATION_TOLERANCE);
    
  } else if (commonScaleCanBeProfiled) {
    profileCommonScale(regression, cache);
  }
}

void profileCommonScale(SEXP regression, MERCache* cache)
{
  const int* dims      = DIMS_SLOT(regression);
  double*    deviances = DEV_SLOT(regression);
    
  double   MLDegreesOfFreedom = cache->priorCommonScaleDegreesOfFreedom + (double) dims[n_POS];
  double REMLDegreesOfFreedom = MLDegreesOfFreedom - (double) dims[p_POS];
  
  double a = cache->totalSumOfSquares + 2.0 * (cache->mTwoExponentialTermConstantPart + cache->mTwoExponentialTermVaryingPart);
  double b, df;
  
  switch (cache->commonScaleOptimization) {
    case CSOT_LINEAR:
      // (sigma^2)^-(df/2) * exp(-0.5 * a / sigma^2)
      // sigma^2_hat = a / df
      deviances[sigmaML_POS]   = sqrt(a /   MLDegreesOfFreedom);
      deviances[sigmaREML_POS] = sqrt(a / REMLDegreesOfFreedom);
      break;
    case CSOT_QUADRATIC_SIGMA:
      // (sigma^2)^-(df/2) * exp(-0.5 * a / sigma^2 - b / sigma)
      // sigma_hat = 0.5 * (b + sqrt(b^2 + 4 * a * df)) / df
      b = cache->mOneExponentialTermConstantPart + cache->mOneExponentialTermVaryingPart;
      df =   MLDegreesOfFreedom;
      deviances[sigmaML_POS]   = 0.5 * (b + sqrt(b * b + 4.0 * a * df)) / df;
      df = REMLDegreesOfFreedom;
      deviances[sigmaREML_POS] = 0.5 * (b + sqrt(b * b + 4.0 * a * df)) / df;
      break;
    case CSOT_QUADRATIC_SIGMA_SQ:
      // (sigma^2)^(-df/2) * exp(-0.5 * a / sigma^2 - b * sigma^2)
      // sigma^2_hat = (sqrt(df^2 + 8 * a * b) - df) / (4 * b)
      b = cache->twoExponentialTermConstantPart + cache->twoExponentialTermVaryingPart;
      df =   MLDegreesOfFreedom;
      deviances[sigmaML_POS]   = sqrt(0.25 * (sqrt(df * df + 8.0 * a * b) - df) / b);
      df = REMLDegreesOfFreedom;
      deviances[sigmaREML_POS] = sqrt(0.25 * (sqrt(df * df + 8.0 * a * b) - df) / b);
      break;
    default: break;
  }
}

static void updateRegressionForNewCommonScale(SEXP regression, MERCache* cache);

double performOneStepOfNewtonsMethodForCommonScale(SEXP regression, MERCache* cache)
{
  double* deviances = DEV_SLOT(regression);
  int*    dims      = DIMS_SLOT(regression);
  
  double firstDerivative, secondDerivative;
  
  getCommonScaleDerivatives(regression, cache, &firstDerivative, &secondDerivative);
  
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

static void updateRegressionForNewCommonScale(SEXP regression, MERCache* cache)
{
  // the update is only required if the total sum of squares term depends on the common scale,
  // itself being a question of whether or not the unmodeled coefficient prior is
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  if (PRIOR_TYPE_SLOT(unmodeledCoefPrior)                         != PRIOR_TYPE_DIRECT ||
      PRIOR_FAMILIES_SLOT(unmodeledCoefPrior)[0]                  != PRIOR_FAMILY_GAUSSIAN ||
      getCommonScaleBit(PRIOR_SCALES_SLOT(unmodeledCoefPrior)[0]) != PRIOR_COMMON_SCALE_FALSE)
    return;
  
  const int* dims      = DIMS_SLOT(regression);
  double*    deviances = DEV_SLOT(regression);
  
  int numUnmodeledCoefs = dims[p_POS];
  
  // we need to refactor (X'X - Rzx'Rzx + sigma^2 / sigma_beta^2 * I)
  double* lowerRightBlockRightFactorization = RX_SLOT(regression);
  
  // recover the cached version of X'X - Rzx'Rzx
  Memcpy(lowerRightBlockRightFactorization, (const double*) cache->downdatedDenseCrossproduct,
         numUnmodeledCoefs * numUnmodeledCoefs);
  
  
  addGaussianContributionToDenseBlock(regression, lowerRightBlockRightFactorization,
                                      deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS]);
  
  int choleskyResult = getDenseCholeskyDecomposition(lowerRightBlockRightFactorization, numUnmodeledCoefs, TRIANGLE_TYPE_UPPER);
  
  if (choleskyResult > 0) error("Leading minor %d of downdated X'X is not positive definite.", choleskyResult);
  if (choleskyResult < 0) error("Illegal argument %d to cholesky decomposition (dpotrf).", -choleskyResult);
  
  deviances[ldRX2_POS] = 0.0;
  for (int j = 0; j < numUnmodeledCoefs; ++j) {
    deviances[ldRX2_POS] += 2.0 * log(lowerRightBlockRightFactorization[j * (numUnmodeledCoefs + 1)]);
  }
  
  // at this point, we have the correct Rx stored. now we need to
  // compute the new half projection and the new sum of squares
  double* unmodeledCoefProjection = cache->unmodeledCoefProjection;
  
  // copy in (X'Y - Rzx' theta half projection); beta half projection is Rx^-1 times that
  Memcpy(unmodeledCoefProjection, (const double*) cache->downdatedDenseResponseRotation,
         numUnmodeledCoefs);
  
  
  int i_one = 1;
  // solve A'x = b for A an Upper triangular, Tranposed, Non-unit matrix
  F77_CALL(dtrsv)("U", "T", "N",
                  &numUnmodeledCoefs,
                  lowerRightBlockRightFactorization,
                  &numUnmodeledCoefs,
                  unmodeledCoefProjection,
                  &i_one);
  
  // now update the sums of squares
  double newSumOfSquares = getSumOfSquares(unmodeledCoefProjection, numUnmodeledCoefs);
  cache->totalSumOfSquares -= newSumOfSquares - cache->unmodeledCoefProjectionSumOfSquares;
  cache->unmodeledCoefProjectionSumOfSquares = newSumOfSquares;
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
//   f'(x) = a * N + b * SS + c * || Rx^-1 beta.tilde ||^2
//   f"(x) = d * N + e * SS + f * || Rx^-1 beta.tilde ||^2 + g * || Rx^-T Rx^-1 beta.tilde ||^2
//
// other terms relating to the other priors might also exist, but are thankfully
// polynomial in the common scale
void getCommonScaleDerivatives(SEXP regression, MERCache* cache, double* firstDerivative, double* secondDerivative)
{
  const int*    dims      = DIMS_SLOT(regression);
  const double* deviances = DEV_SLOT(regression);
  
  double sigma    = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  double sigma_sq = sigma * sigma;
  double sigma_cu = sigma * sigma_sq;
  double sigma_fo = sigma_sq * sigma_sq;
  
  // handle the polynomial term first
  double degreesOfFreedom = cache->priorCommonScaleDegreesOfFreedom;
  degreesOfFreedom += (double) (dims[n_POS] - (dims[isREML_POS] ? dims[p_POS] : 0));
  
  *firstDerivative  = -degreesOfFreedom / sigma;
  *secondDerivative =  degreesOfFreedom / sigma_sq;
  
  // handle the sum of squares and likelihood part directly
  *firstDerivative  +=       cache->totalSumOfSquares / sigma_cu;
  *secondDerivative -= 3.0 * cache->totalSumOfSquares / sigma_fo;
  
  // if the unmodeled coef prior is not on the common scale, further derivatives are required
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  if (PRIOR_TYPE_SLOT(unmodeledCoefPrior) == PRIOR_TYPE_DIRECT &&
      PRIOR_FAMILIES_SLOT(unmodeledCoefPrior)[0] == PRIOR_FAMILY_GAUSSIAN &&
      getCommonScaleBit(PRIOR_SCALES_SLOT(unmodeledCoefPrior)[0]) == PRIOR_COMMON_SCALE_FALSE)
  {
    getDerivativesOfSumOfSquares(regression, cache, firstDerivative, secondDerivative);
  }
  
  
  // now for prior parts
  // exp(-a / sigma^2) => 2a * / sigma^3, -6a / sigma^4
  double a = cache->mTwoExponentialTermConstantPart + cache->mTwoExponentialTermVaryingPart;
  *firstDerivative  += 2.0 * a / sigma_cu;
  *secondDerivative -= 6.0 * a / sigma_fo;
  
  // exp(-a / sigma) => a / sigma^2, -2a / sigma^3
  a = cache->mOneExponentialTermConstantPart + cache->mOneExponentialTermVaryingPart;
  *firstDerivative  +=       a / sigma_sq;
  *secondDerivative -= 2.0 * a / sigma_cu;
  
  // exp(-a * sigma^2) => -2a * sigma, -2a
  a = cache->twoExponentialTermConstantPart + cache->twoExponentialTermVaryingPart;
  *firstDerivative  -= 2.0 * a * sigma;
  *secondDerivative -= 2.0 * a;
  
  // exp(-a * sigma) => -a, 0
  *firstDerivative -= cache->oneExponentialTermConstantPart + cache->oneExponentialTermVaryingPart;
}

static void getDerivativesOfSumOfSquares(SEXP regression, MERCache* cache,
                                         double* firstDerivative, double* secondDerivative)
{
  
  const int* dims      = DIMS_SLOT(regression);
  double*    deviances = DEV_SLOT(regression);
  
  int i_one = 1;
  double d_one = 1.0;
  
  int numUnmodeledCoefs = dims[p_POS];
  
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(unmodeledCoefPrior) + 1; // skip over the log det of the covar, not needed here
  unsigned int numHyperparameters = LENGTH(GET_SLOT(unmodeledCoefPrior, blme_prior_hyperparametersSym)) - 1;
  
  
  // take Rx and get Rx^-1
  const double* lowerRightFactor = RX_SLOT(regression);
  double rightFactorInverse[numUnmodeledCoefs * numUnmodeledCoefs]; // Rx^-1
  invertUpperTriangularMatrix(lowerRightFactor, numUnmodeledCoefs, rightFactorInverse);
  
  // calculate Lbeta^-1 * Rx^-1
  int factorIsTriangular = TRUE;
  if (numHyperparameters == 1) {
    // multiply by a scalar
    // printMatrix(lowerRightFactor, numUnmodeledCoefs, numUnmodeledCoefs);
    for (int col = 0; col < numUnmodeledCoefs; ++col) {
      int offset = col * numUnmodeledCoefs;
      for (int row = 0; row <= col; ++row) {
        rightFactorInverse[offset++] *= hyperparameters[0];
      }
    }
    // printMatrix(rightFactorInverse, numUnmodeledCoefs, numUnmodeledCoefs);
  } else if (numHyperparameters == numUnmodeledCoefs) {
    // left multiply by a diagonal matrix
    const double *diagonal = hyperparameters;
    
    for (int col = 0; col < numUnmodeledCoefs; ++col) {
      int offset = col * numUnmodeledCoefs;
      for (int row = 0; row <= col; ++row) {
        rightFactorInverse[offset++] *= diagonal[row];
      }
    }
  } else {
    const double* priorLeftFactorInverse = hyperparameters;
    // want L * R
    // Left multiply, Lower triangluar matrix, No-transpose, Non-unit
    F77_CALL(dtrmm)("L", "L", "N", "N", &numUnmodeledCoefs, &numUnmodeledCoefs, &d_one,
                    (double*) priorLeftFactorInverse, &numUnmodeledCoefs,
                    rightFactorInverse, &numUnmodeledCoefs);
    factorIsTriangular = FALSE;
  }
  
  double projectionRotation[numUnmodeledCoefs];
  Memcpy(projectionRotation, (const double *) cache->unmodeledCoefProjection, numUnmodeledCoefs);
  
  // this step corresponds to Rx^-1 * unmodeled coef projection
  if (factorIsTriangular) {
    // X := A x, A triangular
    F77_CALL(dtrmv)("Upper triangular", "Non transposed", "Non unit diagonal",
                    &numUnmodeledCoefs, rightFactorInverse, &numUnmodeledCoefs,
                    projectionRotation, &i_one);
  } else {
    applyMatrixToVector(rightFactorInverse, numUnmodeledCoefs, numUnmodeledCoefs, FALSE,
                        projectionRotation, projectionRotation);
  }
  
  double firstRotationSumOfSquares = getSumOfSquares(projectionRotation, numUnmodeledCoefs);
  
  // now for Rx^-T Rx^-1 * modeled coef projection
  if (factorIsTriangular) {
    // X: = A' x, A triangular
    F77_CALL(dtrmv)("Upper triangular", "Transposed", "Non unit diagonal",
                    &numUnmodeledCoefs, rightFactorInverse, &numUnmodeledCoefs,
                    projectionRotation, &i_one);
  } else {
    applyMatrixToVector(rightFactorInverse, numUnmodeledCoefs, numUnmodeledCoefs, TRUE,
                        projectionRotation, projectionRotation);
  }
  
  double secondRotationSumOfSquares = getSumOfSquares(projectionRotation, numUnmodeledCoefs);
  
  
  // in general, DoF depends on unmodeled coefficient prior scale, as we can get back those
  // lost DoF. However, in that case we can't get here, where optimization is required.
  double sigma    = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  double sigma_sq = sigma * sigma;
  
  *firstDerivative  -= firstRotationSumOfSquares / sigma;
  *secondDerivative += 3.0 * firstRotationSumOfSquares / sigma_sq + 4.0 * secondRotationSumOfSquares;
  
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
      // want UU'
      singleTriangularMatrixCrossproduct(rightFactorInverse, numUnmodeledCoefs, TRUE,
                                         TRIANGLE_TYPE_UPPER, crossproduct);
    } else {
      singleMatrixCrossproduct(rightFactorInverse, numUnmodeledCoefs, numUnmodeledCoefs,
                               crossproduct, TRUE, TRIANGLE_TYPE_UPPER);
    }
    double firstOrderTrace  = 0.0;
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
    
    *firstDerivative  -=  sigma * firstOrderTrace;
    *secondDerivative += -firstOrderTrace + 2.0 * sigma_sq * secondOrderTrace;
  }
}
