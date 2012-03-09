#include "unmodeledCoefficientPrior.h"

#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "Syms.h"
#include "lmer.h"
#include "blmer.h"
#include "lmer_common.h"
#include "matrix.h"
#include "multivariateNormal.h"

#include "common_inlines.h"

#include "__lmmMerCache.h" // defines MERCache

double calculateGaussianDeviance(SEXP prior, double commonScale,
                                 const double *unmodeledCoefficients, int numUnmodeledCoefs);
double calculateMVTDeviance(SEXP prior, double commonScale,
                            const double *unmodeledCoefficients, int numUnmodeledCoefs);
double calculateUnmodeledCoefficientDeviance(SEXP prior, double commonScale,
                                             const double *unmodeledCoefficients, int numUnmodeledCoefs)
{
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);
  
  if (priorType == PRIOR_TYPE_NONE) return(0.0);
  
  priorFamily_t priorFamily = PRIOR_FAMILIES_SLOT(prior)[0];
  switch (priorFamily) {
    case PRIOR_FAMILY_GAUSSIAN:
      return(calculateGaussianDeviance(prior, commonScale, unmodeledCoefficients, numUnmodeledCoefs));
      break;
    case PRIOR_FAMILY_MVT:
      return(calculateMVTDeviance(prior, commonScale, unmodeledCoefficients, numUnmodeledCoefs));
      break;
    default:
      break;
  }
  return(0.0);
}

double calculateGaussianDeviance(SEXP prior, double commonScale,
                                 const double *parameters, int numParameters)
{  
  priorScale_t scale  = PRIOR_SCALES_SLOT(prior)[0];
  
  double *hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  int numHyperparameters = LENGTH(GET_SLOT(prior, blme_prior_hyperparametersSym)) - 1;
  
  double logDetCov = *hyperparameters++;
  if (scale == PRIOR_SCALE_COMMON) logDetCov += ((double) numParameters) * log(commonScale);
  
  double result = 0.0;
  if (numHyperparameters == 1) {
    double sdInverse = hyperparameters[0];
    if (scale == PRIOR_SCALE_COMMON) sdInverse /= sqrt(commonScale);
    
    result = -2.0 * dmvn(parameters, numParameters, NULL,
                         sdInverse, logDetCov, TRUE);
  } else if (numHyperparameters == numParameters) {
    double sdsInverse[numParameters];
    for (int i = 0; i < numParameters; ++i) sdsInverse[i] = hyperparameters[i] / sqrt(commonScale);
    
    result = -2.0 * dmvn2(parameters, numParameters, NULL,
                          sdsInverse, logDetCov, TRUE);
  } else if (numHyperparameters == 2 * numParameters * numParameters) {
    int covLength = numParameters * numParameters;
    double covInverse[covLength];
    hyperparameters += covLength; // skip over the left factor
    for (int i = 0; i < covLength; ++i) covInverse[i] = hyperparameters[i] / commonScale;
    
    result = -2.0 * dmvn3(parameters, numParameters, NULL,
                          covInverse, logDetCov, TRUE);
  } else error("Internal error: for a normal prior there are %d hyperparameters but %d coefficients.",
               numHyperparameters + 1, numParameters);
  
  
  return (result);
}

double calculateMVTDeviance(SEXP prior, double commonScale, const double *unmodeledCoefficients, int numUnmodeledCoefs)
{
  error("mvt not yet implemented");
  return (0.0);
}

void addGaussianContributionToDenseBlock(SEXP regression, double* lowerRightBlock, double commonScale)
{
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  priorType_t priorType = PRIOR_TYPE_SLOT(unmodeledCoefPrior);
  
  if (priorType != PRIOR_TYPE_DIRECT) return;
  
  priorFamily_t family = PRIOR_FAMILIES_SLOT(unmodeledCoefPrior)[0];
  
  if (family != PRIOR_FAMILY_GAUSSIAN) return;
  
  
  
  int *dims = DIMS_SLOT(regression);
  int numUnmodeledCoefs = dims[p_POS];
  
  double commonVariance = commonScale * commonScale; //DEV_SLOT(regression)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  // commonVariance *= commonVariance;
  
  priorScale_t scale  = PRIOR_SCALES_SLOT(unmodeledCoefPrior)[0];
  double *hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(unmodeledCoefPrior) + 1; // skip over the log det of the covar, not needed here
  int numHyperparameters = LENGTH(GET_SLOT(unmodeledCoefPrior, blme_prior_hyperparametersSym)) - 1;
  
  if (numHyperparameters == 1) {
    // hyperparameters are log(prior.sd^2) -- skipped, 1 / prior.sd
    double additiveFactor = hyperparameters[0] * hyperparameters[0];
    
    if (scale == PRIOR_SCALE_ABSOLUTE) additiveFactor *= commonVariance;
    
    // add to diagonal
    for (int i = 0; i < numUnmodeledCoefs; ++i) {
      lowerRightBlock[i * (numUnmodeledCoefs + 1)] += additiveFactor;
    }
  } else if (numHyperparameters == numUnmodeledCoefs) {
    // prior covariance is a diagonal matrix, so we store 1 / sqrt of those elements
    
    if (scale == PRIOR_SCALE_ABSOLUTE) {
      for (int i = 0; i < numUnmodeledCoefs; ++i) {
        lowerRightBlock[i * (numUnmodeledCoefs + 1)] += commonVariance * hyperparameters[i] * hyperparameters[i];
      } 
    } else {
      for (int i = 0; i < numUnmodeledCoefs; ++i) {
        lowerRightBlock[i * (numUnmodeledCoefs + 1)] += hyperparameters[i] * hyperparameters[i];
      } 
    }
  } else {
    // prior covariance is an arbitrary matrix. first p^2 components are the left-factor-inverse.
    // second p^2 components are the full inverse
    int covarianceMatrixLength = numUnmodeledCoefs * numUnmodeledCoefs;
    double *covarianceInverse = hyperparameters + covarianceMatrixLength;
    
    if (scale == PRIOR_SCALE_ABSOLUTE) {
      // just need to copy in the upper right block
      int offset;
      for (int col = 0; col < numUnmodeledCoefs; ++col) {
        offset = col * numUnmodeledCoefs;
        for (int row = 0; row <= col; ++row) {
          lowerRightBlock[offset] += commonVariance * covarianceInverse[offset];
          ++offset;
        }
      }
    } else {
      int offset;
      for (int col = 0; col < numUnmodeledCoefs; ++col) {
        offset = col * numUnmodeledCoefs;
        for (int row = 0; row <= col; ++row) {
          lowerRightBlock[offset] += covarianceInverse[offset];
          ++offset;
        }
      }
    }
  }
  
  //Rprintf("after:\n");
  //printMatrix(lowerRightBlock, numUnmodeledCoefs, numUnmodeledCoefs);
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
void getDerivatives(SEXP regression, MERCache *cache,
                    double *firstDerivative, double *secondDerivative)
{
  int    *dims      = DIMS_SLOT(regression);
  double *deviances = DEV_SLOT(regression);
  
  int i_one = 1;
  double d_one = 1.0;
  
  int numUnmodeledCoefs = dims[p_POS];
  
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  double *hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(unmodeledCoefPrior) + 1; // skip over the log det of the covar, not needed here
  unsigned int numHyperparameters = LENGTH(GET_SLOT(unmodeledCoefPrior, blme_prior_hyperparametersSym)) - 1;
  
  
  // take Rx and get Rx^-1
  const double *lowerRightFactor = RX_SLOT(regression);
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
  double degreesOfFreedom = cache->priorDegreesOfFreedom + (double) (dims[n_POS] - (dims[isREML_POS] ? dims[p_POS] : 0));
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
    
    *firstDerivative  -=  currCommonSd * firstOrderTrace;
    *secondDerivative += -firstOrderTrace + 2.0 * currCommonVariance * secondOrderTrace;
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
