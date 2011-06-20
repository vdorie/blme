#include <Rdefines.h>

#include <Rmath.h>               /* gamma functions */
#include <R_ext/Lapack.h>        /* for Lapack and BLAS */

#include "wishart.h"
#include "matrix.h"

double dwish(const double *w, int dim, double degreesOfFreedom,
             const double *scale, int useLog)
{
  double d_dim = (double) dim;
  
  double logDetScale = getLogDeterminantOfPositiveDefiniteMatrix(scale, dim);
  double denominator = ((d_dim * M_LN2 + logDetScale) * degreesOfFreedom +
                         d_dim * (d_dim - 1.0) * M_LN_SQRT_PI) / 2.0;
  
  for (int i = 1; i <= dim; ++i) {
    denominator += lgammafn((degreesOfFreedom + 1.0 - (double) i) / 2.0);
  }
  
  //denominator += (degreesOfFreedom * (double) dim / 2.0) * M_LN2 +
  //               (((double) dim * ((double) dim - 1.0)) / 2.0) * M_LN_SQRT_PI;
  
  //denominator += (degreesOfFreedom / 2.0) * getLogDeterminantOfPositiveDefiniteMatrix(scale, dim);
  
  double numerator = ((degreesOfFreedom - d_dim - 1.0) / 2.0) * getLogDeterminantOfPositiveDefiniteMatrix(w, dim);
  
  
  double matrixProduct[dim * dim];
  double tempMatrix[dim * dim];
  
  Memcpy(tempMatrix, scale, dim * dim);
  Memcpy(matrixProduct, w, dim * dim);
  
  double exponent = 0.0;
  
  // calculate scale * X = w, X = scale^{-1} w
  // wishart uses trace of result
  int lapackResult = solveSymmetricSystem(tempMatrix, dim, matrixProduct, dim);
  
  if (lapackResult < 0) {
    error("error inverting scale for density of wishart, argument %d illegal", -lapackResult);
  } else if (lapackResult > 0) {
    // matrix was not invertible, exponential term should be e^-\infty, so we leave it at 0
  } else {
    double trace = 0.0;
    for (int i = 0; i < dim; ++i) {
      trace += matrixProduct[(dim + 1) * i];
    }
    exponent += -0.5 * trace;
  }
  
  double result = numerator + exponent - denominator;
  
  if (!useLog) return (exp(result));
  else return(result);
}

// wishart with some common calculations pre-made
double dwish2(const double *w, int dim, double degreesOfFreedom,
             double logScaleDeterminant, const double *scaleInverse, int useLog)
{
  double d_dim = (double) dim;
  double denominator = ((d_dim * M_LN2 + logScaleDeterminant) * degreesOfFreedom +
                         d_dim * (d_dim - 1.0) * M_LN_SQRT_PI) / 2.0;
  
  for (int i = 1; i <= dim; ++i) {
    denominator += lgammafn((degreesOfFreedom + 1.0 - (double) i) / 2.0);
  }
  
  //denominator += (degreesOfFreedom * (double) dim / 2.0) * M_LN2 +
  //               (((double) dim * ((double) dim - 1.0)) / 2.0) * M_LN_SQRT_PI;
  
  //denominator += (degreesOfFreedom / 2.0) * logScaleDeterminant;
  
  double numerator = ((degreesOfFreedom - (double) dim - 1.0) / 2.0) * getLogDeterminantOfPositiveDefiniteMatrix(w, dim);
  
  double matrixProduct[dim * dim];
  
  double one =  1.0, zero = 0.0;
  // "L"eft side, such that we multiply by the symmetric matrix on the left (both are symmetric here)
  // "L"ower triangle of symmetric matrix (both triangles are supplied)
  F77_CALL(dsymm)("L", "U", &dim, &dim, &one, scaleInverse, &dim,
           w, &dim, &zero, matrixProduct, &dim);
  
  double exponent = 0.0;
  for (int i = 0; i < dim; ++i) {
    exponent += matrixProduct[(dim + 1) * i];
  }
  exponent *= -0.5;
  
  double result = numerator + exponent - denominator;
  
  if (!useLog) return(exp(result));
  return(result);
}

double diwish(const double *w, int dim, double degreesOfFreedom,
              const double *inverseScale, int useLog)
{
  double denominator = 0.0;
  for (int i = 1; i <= dim; ++i) denominator += lgammafn((degreesOfFreedom + 1.0 - (double) i) / 2.0);
  
  denominator += (degreesOfFreedom * (double) dim / 2.0) * M_LN2 +
                 (((double) dim * ((double) dim - 1.0)) / 2.0) * M_LN_SQRT_PI;
  
  denominator += ((degreesOfFreedom + (double) dim + 1.0) / 2.0) * getLogDeterminantOfPositiveDefiniteMatrix(w, dim);
  
  double numerator = (degreesOfFreedom / 2.0) * getLogDeterminantOfPositiveDefiniteMatrix(inverseScale, dim);
  
  double tempMatrix[dim * dim];
  double matrixProduct[dim * dim];
  
  Memcpy(tempMatrix, w, dim * dim);
  Memcpy(matrixProduct, inverseScale, dim * dim);
  
  double exponent = 0.0;
  
  // calculate w * X = inverseScale, X = w^{-1} inverseScale
  // wishart uses trace of result
  int lapackResult = solveSymmetricSystem(tempMatrix, dim, matrixProduct, dim);
  
  if (lapackResult < 0) {
    error("error inverting parameter matrix for density of inverse wishart, argument %d illegal", -lapackResult);
  } else if (lapackResult > 0) {
    // matrix was not invertible, exponential term should be e^-\infty, so we leave it at 0
  } else {
    double trace = 0.0;
    for (int i = 0; i < dim; ++i) {
      trace += matrixProduct[(dim + 1) * i];;
    }
    
    exponent += -0.5 * trace;
  }
  
  double result = numerator + exponent - denominator;
  
  if (!useLog) return (exp(result));
  else return(result);
}

double diwish2(const double *wInverse, int dim, double degreesOfFreedom,
               double logInverseScaleDeterminant, const double *inverseScale, int useLog)
{
  double denominator = 0.0;
  for (int i = 1; i <= dim; ++i) denominator += lgammafn((degreesOfFreedom + 1.0 - (double) i) / 2.0);
  
  denominator += (degreesOfFreedom * (double) dim / 2.0) * M_LN2 +
                 (((double) dim * ((double) dim - 1.0)) / 2.0) * M_LN_SQRT_PI;
  
  denominator += -((degreesOfFreedom + (double) dim + 1.0) / 2.0) * getLogDeterminantOfPositiveDefiniteMatrix(wInverse, dim);
  
  double numerator = (degreesOfFreedom / 2.0) * logInverseScaleDeterminant;
  

  
  double matrixProduct[dim * dim];
  
  double one =  1.0, zero = 0.0;
  // "L"eft side, such that we multiply by the symmetric matrix on the left (both are symmetric here)
  // "L"ower triangle of symmetric matrix (both triangles are supplied)
  F77_CALL(dsymm)("L", "U", &dim, &dim, &one, inverseScale, &dim,
           wInverse, &dim, &zero, matrixProduct, &dim);
  
  double exponent = 0.0;
  for (int i = 0; i < dim; ++i) {
    exponent += matrixProduct[(dim + 1) * i];;
  }
    
  exponent *= -0.5;
  
  double result = numerator + exponent - denominator;
  
  if (!useLog) return (exp(result));
  else return(result);
}
