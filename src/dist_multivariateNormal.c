#include <Rdefines.h>

#include <Rmath.h>               /* gamma functions */
#include <R_ext/Lapack.h>        /* for Lapack and BLAS */

#include "dist_multivariateNormal.h"
#include "matrix.h"

/**
 * Density of a multivariate normal with given mean and covariance
 *
 * @param x            normally distributed vector for which the density is desired
 * @param dim          length of x, the mean, and the sides of the covar
 * @param mean         mean of x. can be NULL
 * @param sdInverse    1 / single standard deviation to be used for all (iid)
 * @param sdsInverse   1 / vector of sds (indep)
 * @param covInverse   arbitrary covariance matrix
 * @param logDetCov    the log of the determinant of the covariance matrix
 * @param useLog       whether or not the result is returned as the logarithm of the density
 *
 * @return the density or log density of x
 */
double dmvn(const double *x, int dim, const double *mean, double sdInverse,
            double logDetCov, int useLog)
{
  double exponent = 0.0;
  if (mean != NULL) {
    for (int i = 0; i < dim; ++i) {
      double x_i = x[i] - mean[i];
      exponent += x_i * x_i;
    }
  } else {
    for (int i = 0; i < dim; ++i) exponent += x[i] * x[i];
  }
  exponent *= -0.5 * sdInverse * sdInverse;
  
  double result = exponent - 0.5 * logDetCov - ((double) dim) * (0.5 * M_LN2 + M_LN_SQRT_PI);
  
  if (!useLog) return(exp(result));
  return(result);
}
double dmvn2(const double *x, int dim, const double *mean, const double* sdsInverse,
             double logDetCov, int useLog)
{
  double exponent = 0.0;
  if (mean != NULL) {
    for (int i = 0; i < dim; ++i) {
      double x_i = x[i] - mean[i];
      exponent += x_i * x_i;
    }
  } else {
    for (int i = 0; i < dim; ++i) exponent += x[i] * x[i] * sdsInverse[i] * sdsInverse[i];
  }
  exponent *= -0.5;
  
  double result = exponent - 0.5 * logDetCov - ((double) dim) * (0.5 * M_LN2 + M_LN_SQRT_PI);
  
  if (!useLog) return(exp(result));
  return(result);
}
double dmvn3(const double *x, int dim, const double *mean, const double* covInverse,
             double logDetCov, int useLog)
{
  double centeredX[dim];
  double temp[dim];
  
  if (mean != NULL) {
    for (int i = 0; i < dim; ++i) centeredX[i] = x[i] - mean[i];
    x = centeredX;
  }
  
  applyMatrixToVector(covInverse, dim, dim, TRUE, x, temp);
  
  double exponent = 0.0;
  for (int i = 0; i < dim; ++i) exponent += x[i] * temp[i];
  exponent *= -0.5;
  
  double result = exponent - 0.5 * logDetCov - ((double) dim) * (0.5 * M_LN2 + M_LN_SQRT_PI);
  
  if (!useLog) return(exp(result));
  return(result);
}
