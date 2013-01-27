#include "unmodeledCoefficientPrior.h"

#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "Syms.h"
#include "lmer.h"
#include "blmer.h"
#include "lmer_common.h"
#include "matrix.h"

#include "common_inlines.h"

#include "__lmmMerCache.h" // defines MERCache

double getGaussianDevianceVaryingPart(SEXP prior, double commonScale,
                                      const double *unmodeledCoefficients, int numUnmodeledCoefs);

double getMVTDevianceVaryingPart(SEXP prior, double commonScale,
                                 const double *unmodeledCoefficients, int numUnmodeledCoefs);

double getGaussianDevianceConstantPart(SEXP prior, int numUnmodeledCoefs);
double getMVTDevianceConstantPart(SEXP prior, int numUnmodeledCoefs);

double getUnmodeledCoefficientDevianceVaryingPart(SEXP prior, double commonScale,
                                                  const double *unmodeledCoefficients, int numUnmodeledCoefs)
{
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);
  
  if (priorType != PRIOR_TYPE_DIRECT) return(0.0);
  
  priorFamily_t priorFamily = PRIOR_FAMILIES_SLOT(prior)[0];
  switch (priorFamily) {
    case PRIOR_FAMILY_GAUSSIAN:
      return(getGaussianDevianceVaryingPart(prior, commonScale, unmodeledCoefficients, numUnmodeledCoefs));
      break;
    case PRIOR_FAMILY_MVT:
      return(getMVTDevianceVaryingPart(prior, commonScale, unmodeledCoefficients, numUnmodeledCoefs));
      break;
    default:
      break;
  }
  return(0.0);
}

double getUnmodeledCoefficientDevianceConstantPart(SEXP prior, int numUnmodeledCoefs)
{
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);
  
  if (priorType != PRIOR_TYPE_DIRECT) return(0.0);
  
  priorFamily_t priorFamily = PRIOR_FAMILIES_SLOT(prior)[0];
  switch (priorFamily) {
    case PRIOR_FAMILY_GAUSSIAN:
      return(getGaussianDevianceConstantPart(prior, numUnmodeledCoefs));
      break;
    case PRIOR_FAMILY_MVT:
      return(getMVTDevianceConstantPart(prior, numUnmodeledCoefs));
      break;
    default:
      break;
  }
  return(0.0);
}

double getUnmodeledCoefficientDensityExponentialPart(SEXP prior, const double* unmodeledCoefficients, int numUnmodeledCoefs)
{
  if (PRIOR_TYPE_SLOT(prior) != PRIOR_TYPE_DIRECT ||
      PRIOR_FAMILIES_SLOT(prior)[0] != PRIOR_FAMILY_GAUSSIAN) return(0.0);
  
  const double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior) + 1;
  int numHyperparameters  = LENGTH(GET_SLOT(prior, blme_prior_hyperparametersSym)) - 1; // pop off one for the log determinant of the covariance
  
  if (numHyperparameters == 1) {
    double varInverse = hyperparameters[0] * hyperparameters[0];
    
    return(getSumOfSquares(unmodeledCoefficients, numUnmodeledCoefs) * varInverse);
  } else if (numHyperparameters == numUnmodeledCoefs) {
    const double* sdsInverse = hyperparameters;
    
    double result = 0.0;
    for (int i = 0; i < numUnmodeledCoefs; ++i) {
      result += unmodeledCoefficients[i] * unmodeledCoefficients[i] * sdsInverse[i] * sdsInverse[i];
    }
    
    return(result);
  } else if (numHyperparameters == 2 * numUnmodeledCoefs * numUnmodeledCoefs) {
    int covLength = numUnmodeledCoefs * numUnmodeledCoefs;
    // hyperparameters contains: ( logDetCov, leftFactor, covInv ). skip left factor
    const double* covInverse = hyperparameters + covLength;

    double tempVector[numUnmodeledCoefs];
    
    applyMatrixToVector(covInverse, numUnmodeledCoefs, numUnmodeledCoefs, TRUE, unmodeledCoefficients, tempVector);
    
    double result = 0.0;
    for (int i = 0; i < numUnmodeledCoefs; ++i) {
      result += unmodeledCoefficients[i] * tempVector[i];
    }
    
    return(result);
  }
  
  return(0.0);
}
  

double getGaussianDevianceVaryingPart(SEXP prior, double commonScale,
                                      const double *parameters, int numParameters)
{  
  priorCommonScale_t useCommonScale = getCommonScaleBit(PRIOR_SCALES_SLOT(prior)[0]);
  
  
  double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior) + 1;
  int numHyperparameters  = LENGTH(GET_SLOT(prior, blme_prior_hyperparametersSym)) - 1; // pop off one for the log determinant of the covariance
  
  double result = 0.0;
  if (useCommonScale) result = ((double) numParameters) * log(commonScale);
  
  if (numHyperparameters == 1) {
    double varInverse = hyperparameters[0] * hyperparameters[0];
    if (useCommonScale) varInverse /= commonScale;
    
    result += getSumOfSquares(parameters, numParameters) * varInverse;
    
  } else if (numHyperparameters == numParameters) {
    double* sdsInverse;
    double temp[numParameters];
    if (useCommonScale) {
      for (int i = 0; i < numParameters; ++i) temp[i] = hyperparameters[i] / sqrt(commonScale);
      sdsInverse = temp;
    } else {
      sdsInverse = hyperparameters;
    }
    
    for (int i = 0; i < numParameters; ++i) {
      result += parameters[i] * parameters[i] * sdsInverse[i] * sdsInverse[i];
    }
    
  } else if (numHyperparameters == 2 * numParameters * numParameters) {
    int covLength = numParameters * numParameters;
    double* covInverse;
    double tempMatrix[covLength];
    double tempVector[numParameters];
    
    hyperparameters += covLength; // hyperparameters contains: ( logDetCov, leftFactor, covInv ). skip left factor
    if (useCommonScale) {
      for (int i = 0; i < covLength; ++i) tempMatrix[i] = hyperparameters[i] / commonScale;
      covInverse = tempMatrix;
    } else {
      covInverse = hyperparameters;
    }
    
    
    applyMatrixToVector(covInverse, numParameters, numParameters, TRUE, parameters, tempVector);
    
    for (int i = 0; i < numParameters; ++i) {
      result += parameters[i] * tempVector[i];
    }
    
  } else error("Internal error: for a normal prior there are %d hyperparameters but %d coefficients.",
               numHyperparameters + 1, numParameters);
  
  
  return (result);
}

double getGaussianDevianceConstantPart(SEXP prior, int numParameters)
{
  double* hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  double logDetCov = *hyperparameters++;
  
  return (logDetCov + ((double) numParameters) * (M_LN2 + 2.0 * M_LN_SQRT_PI));
}

double getMVTDevianceVaryingPart(SEXP prior, double commonScale, const double *unmodeledCoefficients, int numUnmodeledCoefs)
{
  error("mvt not yet implemented");
  return (0.0);
}

double getMVTDevianceConstantPart(SEXP prior, int numUnmodeledCoefs)
{
  error("mvt not yet implemented");
  return (0.0);
}

double getUnmodeledCoefficientPriorCommonScaleDegreesOfFreedom(SEXP prior, int numUnmodeledCoefs)
{
  if (PRIOR_TYPE_SLOT(prior) == PRIOR_TYPE_DIRECT &&
      PRIOR_FAMILIES_SLOT(prior)[0] == PRIOR_FAMILY_GAUSSIAN &&
      getCommonScaleBit(PRIOR_SCALES_SLOT(prior)[0]) == PRIOR_COMMON_SCALE_TRUE)
  {
    return((double) numUnmodeledCoefs);
  }
  
  return(0.0);
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
  
  priorCommonScale_t useCommonScale = getCommonScaleBit(PRIOR_SCALES_SLOT(unmodeledCoefPrior)[0]);
  
  double *hyperparameters = PRIOR_HYPERPARAMETERS_SLOT(unmodeledCoefPrior) + 1; // skip over the log det of the covar, not needed here
  int numHyperparameters = LENGTH(GET_SLOT(unmodeledCoefPrior, blme_prior_hyperparametersSym)) - 1;
  
  if (numHyperparameters == 1) {
    // hyperparameters are log(prior.sd^2) -- skipped, 1 / prior.sd
    double additiveFactor = hyperparameters[0] * hyperparameters[0];
    
    if (!useCommonScale) additiveFactor *= commonVariance;
    
    // add to diagonal
    for (int i = 0; i < numUnmodeledCoefs; ++i) {
      lowerRightBlock[i * (numUnmodeledCoefs + 1)] += additiveFactor;
    }
  } else if (numHyperparameters == numUnmodeledCoefs) {
    // prior covariance is a diagonal matrix, so we store 1 / sqrt of those elements
    
    if (!useCommonScale) {
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
    
    if (!useCommonScale) {
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
