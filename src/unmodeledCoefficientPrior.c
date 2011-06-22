#include "unmodeledCoefficientPrior.h"
#include "blmer.h"
#include "multivariateNormal.h"
#include "lmer.h"

#include <Rmath.h>
#include "Syms.h"

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
  } else if (numHyperparameters == numParameters * numParameters) {
    int covLength = numParameters * numParameters;
    double covInverse[covLength];
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

void addGaussianContributionToDenseBlock(SEXP regression, double *lowerRightBlock)
{
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  priorType_t priorType = PRIOR_TYPE_SLOT(unmodeledCoefPrior);
  
  if (priorType != PRIOR_TYPE_DIRECT) return;
  
  priorFamily_t family = PRIOR_FAMILIES_SLOT(unmodeledCoefPrior)[0];
  
  if (family != PRIOR_FAMILY_GAUSSIAN) return;
  
  
  
  int *dims = DIMS_SLOT(regression);
  int numUnmodeledCoefs = dims[p_POS];
  
  double commonVariance = DEV_SLOT(regression)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  commonVariance *= commonVariance;
  
  priorScale_t scale  = PRIOR_SCALES_SLOT(unmodeledCoefPrior)[0];
  double *hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(unmodeledCoefPrior) + 1; // skip over the log det of the covar, not needed here
  int numHyperparameters = LENGTH(GET_SLOT(unmodeledCoefPrior, blme_prior_hyperparametersSym)) - 1;
  
  if (numHyperparameters == 1) {
    // hyperparameters are log(prior.sd^2), 1 / prior.sd
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
