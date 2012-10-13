#include "covariancePrior.h"

#include <Rmath.h>
#include "Syms.h"

#include "blmer.h"
#include "parameters.h"
#include "util.h"
#include "wishart.h"

// forward declarations
static double calculateCorrelationDeviance(SEXP prior, double commonScale, const double *parameters,
                                           int factorDimension);
static double calculateSpectralDeviance(SEXP prior, double commonScale, const double *parameters,
                                        int factorDimension);
static double calculateDirectDeviance(SEXP prior, double commonScale, const double *parameters,
                                      int factorDimension);
static void setCorrelationConstraints(double *constraints, int dim);
static void setSpectralConstraints(double *constraints, int dim);
static void setSTConstraints(double *constraints, int dim);

// exported functions
double calculateCovarianceDeviance(SEXP prior, double commonScale, const double *parameters,
                                   int factorDimension)
{
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);
  
  double result = 0.0;
  commonScale = 1.0;
  
  switch (priorType) {
    case PRIOR_TYPE_CORRELATION:
      result = calculateCorrelationDeviance(prior, commonScale, parameters, factorDimension);
      break;
    case PRIOR_TYPE_SPECTRAL:
      result = calculateSpectralDeviance(prior, commonScale, parameters, factorDimension);
      break;
    case PRIOR_TYPE_DIRECT:
      result = calculateDirectDeviance(prior, commonScale, parameters, factorDimension);
    default:
      break;
  }
  
  return(result);
}

void setCovarianceConstraints(SEXP prior, double *boxConstraints, int factorDimension)
{
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);
  switch (priorType) {
    case PRIOR_TYPE_CORRELATION:
      setCorrelationConstraints(boxConstraints, factorDimension);
      break;
      
    case PRIOR_TYPE_SPECTRAL:
      setSpectralConstraints(boxConstraints, factorDimension);
      break;
    
    case PRIOR_TYPE_DIRECT:
      // direct priors are done using the ST decomp
      // setDirectConstraints(boxConstraints, factorDimension);
      setSTConstraints(boxConstraints, factorDimension);
      break;
    
    default:
      setSTConstraints(boxConstraints, factorDimension);
      break;
      
  }
}


// internal functions

// deviances
static double calculateCorrelationDeviance(SEXP prior, double commonScale, const double *parameters,
                                           int factorDimension)
{
  double deviance = 0.0;
  
  commonScale = sqrt(commonScale);
  
  priorScale_t *scales    = PRIOR_SCALES_SLOT(prior);
  priorFamily_t *familes  = PRIOR_FAMILIES_SLOT(prior);
  double *hyperParameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  for (int i = 0; i < factorDimension; ++i) {
    switch (familes[i]) {
      case (PRIOR_FAMILY_GAMMA):
        {
          double shape = *hyperParameters++;
          double scale = 1 / *hyperParameters++; // we store default as rate, which is the R default. C uses scale
          
          double parameter = *parameters++ * sqrt(commonScale);
          if (*scales++ == PRIOR_SCALE_VARIANCE) parameter *= parameter;
          
          if (parameter == 0.0) {
            if (shape > 1.0) {  // for shape > 1, gamma densities go to 0 at 0, otherwise go to infinity
              return ( INFINITY); // -2 * log(0) = inf
            } else {
              return (-INFINITY); // -2 * log(inf) = -inf
            }
          }
        
          deviance += -2.0 * dgamma(parameter, shape, scale, TRUE);
        }
        
        break;
        
      case (PRIOR_FAMILY_INVGAMMA):
        {
          double shape = *hyperParameters++;
          double scale = *hyperParameters++;
          
          double parameter = *parameters++ * sqrt(commonScale);
          if (*scales++ == PRIOR_SCALE_VARIANCE) parameter *= parameter;
          
          if (parameter == 0.0) return (INFINITY); // -2 * log(0) = inf
        
          deviance += -2.0 * (dgamma(1 / parameter, shape, 1 / scale, TRUE) - 2.0 * log(parameter));
        }
        
        break;
      default:
        break;
    }
  }
  
  switch (familes[factorDimension]) {
    case PRIOR_FAMILY_WISHART:
      {
        double degreesOfFreedom = *hyperParameters++;
        double logScaleDeterminant = *hyperParameters++;
        double *scaleInverse = hyperParameters;
        hyperParameters += factorDimension * factorDimension;
        
        double *centerMatrix = Alloca(factorDimension * factorDimension, double);
        R_CheckStack();
        convertSTToCovarianceWithScale(parameters, factorDimension, centerMatrix,
                                       commonScale);
        
        deviance += -2.0 * dwish2(centerMatrix, factorDimension, degreesOfFreedom, 
                                  logScaleDeterminant, scaleInverse, TRUE);
      }
      break;
    case PRIOR_FAMILY_INVWISHART:
      {
        double degreesOfFreedom = *hyperParameters++;
        double logInverseScaleDeterminant = *hyperParameters++;
        double *inverseScale = hyperParameters;
        hyperParameters += factorDimension * factorDimension;
        
        double *centerMatrixInverse = Alloca(factorDimension * factorDimension, double);
        R_CheckStack();
        convertSTToCovarianceInverseWithScale(parameters, factorDimension, centerMatrixInverse,
                                              1.0 / commonScale);
        
        deviance += -2.0 * diwish2(centerMatrixInverse, factorDimension, degreesOfFreedom, 
                                   logInverseScaleDeterminant, inverseScale, TRUE);
      }
      break;
    default:
      break;
  }
    
  parameters += factorDimension * (factorDimension + 1) / 2;
  
  return (deviance);
}

static double calculateSpectralDeviance(SEXP prior, double commonScale, const double *parameters, int factorDimension)
{
  double deviance = 0.0;
  
  priorScale_t *scales = PRIOR_SCALES_SLOT(prior);
  priorFamily_t *familes = PRIOR_FAMILIES_SLOT(prior);
  double *hyperParameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  for (int i = 0; i < factorDimension; ++i) {
    switch (familes[i]) {
      case (PRIOR_FAMILY_GAMMA):
        {
          double shape = *hyperParameters++;
          double scale = 1 / *hyperParameters++;
          
          // the raw parameter is an eigen value, which is on the scale of a variance
          double parameter = *parameters++ * commonScale;
          if (*scales++ == PRIOR_SCALE_SD) parameter = sqrt(parameter);
          
          if (parameter == 0.0) {
            if (shape > 1.0) {  // for shape > 1, gamma densities go to 0 at 0, otherwise go to infinity
              return ( INFINITY); // -2 * log(0) = inf
            } else {
              return (-INFINITY); // -2 * log(inf) = -inf
            }
          }
          
          // Rprintf("param = %f, gamma deviance = %f, shape = %f, rate = %f\n", parameter, -2.0 * dgamma(parameter, shape, rate, TRUE), shape, rate);
          
          deviance += -2.0 * dgamma(parameter, shape, scale, TRUE);
        }
        
        break;
      case (PRIOR_FAMILY_INVGAMMA):
        {
          double shape = *hyperParameters++;
          double scale = *hyperParameters++;
          
          // same as above
          double parameter = *parameters++ * commonScale;
          if (*scales++ == PRIOR_SCALE_SD) parameter = sqrt(parameter);
          
          if (parameter == 0.0) return (INFINITY); // -2 * log(0) = inf
          
          // Rprintf("param = %f, inv-gamma deviance = %f, shape = %f, scale = %f\n", parameter, -2.0 * (dgamma(1 / parameter, shape, scale, TRUE) - 2.0 * log(parameter)), shape, scale);
          
          deviance += -2.0 * (dgamma(1 / parameter, shape, 1 / scale, TRUE) - 2.0 * log(parameter));
        }
        
        break;

      default:
        break;
    }
  }
  
//  Rprintf("spectral deviance: %f\n", deviance);
  
  return (deviance);
}

static double calculateDirectDeviance(SEXP prior, double commonScale, const double *parameters,
                                      int factorDimension)
{
  double deviance = 0.0;
  
  priorScale_t *scales    = PRIOR_SCALES_SLOT(prior);
  priorFamily_t *familes  = PRIOR_FAMILIES_SLOT(prior);
  double *hyperParameters = PRIOR_HYPERPARAMETERS_SLOT(prior);
  
  commonScale = sqrt(commonScale);
  
  switch (familes[0]) {
    case (PRIOR_FAMILY_GAMMA):
      {
        double shape = *hyperParameters++;
        double scale = 1 / *hyperParameters++;
          
        // the raw parameter is a standard deviation
        double parameter = *parameters++ * commonScale;
        if (*scales++ == PRIOR_SCALE_VARIANCE) parameter *= parameter;
        
        if (parameter == 0.0) {
          if (shape > 1.0) {  // for shape > 1, gamma densities go to 0 at 0, otherwise go to infinity
            return  (INFINITY); // -2 * log(0) = inf
          } else {
            return (-INFINITY); // -2 * log(inf) = -inf
          }
        }
        
        deviance += -2.0 * dgamma(parameter, shape, scale, TRUE);
      }
      
      break;
    case (PRIOR_FAMILY_INVGAMMA):
      {
        double shape = *hyperParameters++;
        double scale = *hyperParameters++;
          
        // same as above
        double parameter = *parameters++ * commonScale;
        if (*scales++ == PRIOR_SCALE_VARIANCE) parameter *= parameter;
          
        if (parameter == 0.0) return (INFINITY); // -2 * log(0) = inf
                    
        deviance += -2.0 * (dgamma(1.0 / parameter, shape, 1.0 / scale, TRUE) - 2.0 * log(parameter));
      }
      break;
      
    case PRIOR_FAMILY_WISHART:
      {
        double degreesOfFreedom = *hyperParameters++;
        double logScaleDeterminant = *hyperParameters++;
        double *scaleInverse = hyperParameters;
        hyperParameters += factorDimension * factorDimension;
        
        double *covariance = Alloca(factorDimension * factorDimension, double);
        R_CheckStack();
        convertSTToCovarianceWithScale(parameters, factorDimension, covariance,
                                       commonScale * commonScale);
        
        deviance += -2.0 * dwish2(covariance, factorDimension, degreesOfFreedom, 
                                  logScaleDeterminant, scaleInverse, TRUE);
      }
      break;
      
    case PRIOR_FAMILY_INVWISHART:
      {
        double degreesOfFreedom = *hyperParameters++;
        double logInverseScaleDeterminant = *hyperParameters++;
        double *inverseScale = hyperParameters;
        hyperParameters += factorDimension * factorDimension;
        
        double *covarianceInverse = Alloca(factorDimension * factorDimension, double);
        R_CheckStack();
        convertSTToCovarianceInverseWithScale(parameters, factorDimension, covarianceInverse,
                                              1.0 / (commonScale * commonScale));
        
        deviance += -2.0 * diwish2(covarianceInverse, factorDimension, degreesOfFreedom, 
                                   logInverseScaleDeterminant, inverseScale, TRUE);
      }
      break;
      
    default:
      break;
  }
  
  return (deviance);
}

// constraints
static void setCorrelationConstraints(double *constraints, int dim)
{
  // correlation is scales, diagonal, lower triangle.
  // in theory, we should keep the middle matrix positive definite,
  // but we let our prior handle that
        
  // scales and diagonals
  int numParameters = 2 * dim;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = 0.0;
    *constraints++ = R_PosInf;
  }
        
  // lower triangle of positive definite matrix
  numParameters = dim * (dim - 1) / 2;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = R_NegInf;
    *constraints++ = R_PosInf;
  }
}


static void setSpectralConstraints(double *constraints, int dim)
{
  // spectrals have eigenvalues then angles. The angles should be (-pi, pi],
  // but it is doubtful that the optimization algorithm works on the circle so
  // we extend it to the real line with the understanding that any mode is good
  // enough.
  
  // eigenvalues
  int numParameters = dim;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = 0.0;
    *constraints++ = R_PosInf;
  }
  
  // cayley elements
  numParameters = dim * (dim - 1) / 2;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = R_NegInf;
    *constraints++ = R_PosInf;
  }
}

static void setSTConstraints(double *constraints, int dim)
{
  // default is on the cholesky scale, S first then T
  // S is positive, T is unconstrained
  
  // S
  int numParameters = dim;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = 0.0;
    *constraints++ = R_PosInf;
  }
  
  // T
  numParameters = dim * (dim - 1) / 2;
  for (int j = 0; j < numParameters; ++j) {
    *constraints++ = R_NegInf;
    *constraints++ = R_PosInf;
  }
}
