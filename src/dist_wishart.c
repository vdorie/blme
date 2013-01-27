#include "dist_wishart.h"

#include <Rdefines.h>
#include <Rmath.h>               /* gamma functions and finite check */
#include <R_ext/Lapack.h>        /* for Lapack and BLAS */

#include "parameters.h"
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
  // I hope this never happens, but it is possible to use a degenerate scale with a psuedo
  // inverse
  double denominator = (d_dim * M_LN2 * degreesOfFreedom +
                         d_dim * (d_dim - 1.0) * M_LN_SQRT_PI) / 2.0;
  if (R_FINITE(logScaleDeterminant)) denominator += 0.5 * logScaleDeterminant * degreesOfFreedom;
  
  for (int i = 1; i <= dim; ++i) {
    denominator += lgammafn((degreesOfFreedom + 1.0 - (double) i) / 2.0);
  }
  
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
  
  double logInverseScaleDeterminant = getLogDeterminantOfPositiveDefiniteMatrix(inverseScale, dim);
  double numerator = (R_FINITE(logInverseScaleDeterminant) ? 0.5 * degreesOfFreedom * logInverseScaleDeterminant : 0.0);
  
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
  
  double numerator = (R_FINITE(logInverseScaleDeterminant) ? 0.5 * degreesOfFreedom * logInverseScaleDeterminant : 0.0);
  

  
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

double wishartDevianceVaryingPart(const double* st, int dim, double degreesOfFreedom, const double* scaleInverse)
{
  // times 2.0 as we taking the det of only a square root of the desired matrix
  // also not times 2.0 as we might expect for a deviance as there determinant is raised to (stuff / 2)
  double logDet = getLogDeterminantOfSTMatrix(st, dim, NULL);
  if (ISNAN(logDet)) return(R_NaN);
  if (!R_FINITE(logDet)) return(logDet);
  
  double result = -(degreesOfFreedom - (double) dim - 1.0) * 2.0 * logDet;

  int matrixLength = dim * dim;
  double covariance[matrixLength];
  convertSTToCovariance(st, dim, covariance);
  
  // since we have tr(AB') = sum over hadamard product, and as both matrices are symmetric,
  // the trace is given by the following loop
  for (int i = 0; i < matrixLength; ++i) result += scaleInverse[i] * covariance[i];
  
  return(result);
}

double wishartDevianceConstantPart(int dim, double degreesOfFreedom, double logScaleDeterminant)
{
  if (!R_FINITE(logScaleDeterminant)) return(0.0); // for the improper prior
  
  double d_dim = (double) dim;
  
  double result = degreesOfFreedom * (d_dim * M_LN2 + logScaleDeterminant);
  
  result += d_dim * (d_dim - 1.0) * M_LN_SQRT_PI;
  for (int i = 1; i <= dim; ++i) result += 2.0 * lgammafn((degreesOfFreedom + 1.0 - (double) i) / 2.0);
  
  return(result);
}

double wishartDevianceVaryingPartWithScale(const double* st, int dim, double degreesOfFreedom, const double* scaleInverse, double scale)
{
  double logDet = getLogDeterminantOfSTMatrix(st, dim, NULL);
  if (ISNAN(logDet)) return(R_NaN);
  if (!R_FINITE(logDet)) return(logDet);
  
  double d_dim = (double) dim;
  
  // times 2.0 as we taking the det of only a square root of the desired matrix
  // also not times 2.0 as we might expect for a deviance as there determinant is raised to (stuff / 2)
  double result = -(degreesOfFreedom - d_dim - 1.0) * (2.0 * logDet + d_dim * log(scale));
  
  int matrixLength = dim * dim;
  double covariance[matrixLength];
  convertSTToCovariance(st, dim, covariance);
  
  // since we have tr(AB') = sum over hadamard product, and as both matrices are symmetric,
  // the trace is given by the following loop
  double exponential = 0.0;
  for (int i = 0; i < matrixLength; ++i) exponential += scaleInverse[i] * covariance[i];
  
  result += exponential / scale;
  
  return(result);
}



double wishartDevianceVaryingDeterminantPart(const double* st, int dim, double degreesOfFreedom)
{
  return( -(degreesOfFreedom - (double) dim - 1.0) * 2.0 * getLogDeterminantOfSTMatrix(st, dim, NULL));
}

double wishartDevianceVaryingExponentialPart(const double* st, int dim, const double* scaleInverse)
{
  double result = 0.0;
  
  int matrixLength = dim * dim;
  double covariance[matrixLength];
  convertSTToCovariance(st, dim, covariance);
  
  for (int i = 0; i < matrixLength; ++i) result += scaleInverse[i] * covariance[i];
  
  return(result);
}

double inverseWishartDevianceVaryingPart(const double* st, int dim, double degreesOfFreedom, const double* inverseScale)
{
  double logDet = getLogDeterminantOfSTMatrix(st, dim, NULL);
  if (ISNAN(logDet)) return(R_NaN);
  if (!R_FINITE(logDet)) return(logDet);
  
  double result = (degreesOfFreedom + (double) dim + 1.0) * 2.0 * logDet;
  
  int matrixLength = dim * dim;
  double covarianceInverse[matrixLength];
  convertSTToCovarianceInverse(st, dim, covarianceInverse);
  
  for (int i = 0; i < matrixLength; ++i) result += inverseScale[i] * covarianceInverse[i];
  
  return(result);
}

double inverseWishartDevianceVaryingPartWithScale(const double* st, int dim, double degreesOfFreedom, const double* inverseScale, double scale)
{
  double logDet = getLogDeterminantOfSTMatrix(st, dim, NULL);
  if (ISNAN(logDet)) return(R_NaN);
  if (!R_FINITE(logDet)) return(logDet);
  
  double d_dim = (double) dim;
  
  double result = (degreesOfFreedom + d_dim + 1.0) * (2.0 * logDet + d_dim * log(scale));
  
  int matrixLength = dim * dim;
  double covarianceInverse[matrixLength];
  convertSTToCovarianceInverse(st, dim, covarianceInverse);
  
  double exponential = 0.0;
  for (int i = 0; i < matrixLength; ++i) exponential += inverseScale[i] * covarianceInverse[i];
  
  result += exponential / scale;
  
  return(result);
}


double inverseWishartDevianceVaryingDeterminantPart(const double* st, int dim, double degreesOfFreedom)
{
  return((degreesOfFreedom + (double) dim + 1.0) * 2.0 * getLogDeterminantOfSTMatrix(st, dim, NULL));
}

double inverseWishartDevianceVaryingExponentialPart(const double* st, int dim, const double* inverseScale)
{
  double result = 0.0;
  
  int matrixLength = dim * dim;
  double covarianceInverse[matrixLength];
  convertSTToCovarianceInverse(st, dim, covarianceInverse);
  
  for (int i = 0; i < matrixLength; ++i) result += inverseScale[i] * covarianceInverse[i];
  
  return(result);
}

double inverseWishartDevianceConstantPart(int dim, double degreesOfFreedom, double logInverseScaleDeterminant)
{
  if (!R_FINITE(logInverseScaleDeterminant)) return(0.0); // for the improper prior
  
  double d_dim = (double) dim;
  
  double result = degreesOfFreedom * (d_dim * M_LN2 - logInverseScaleDeterminant);
  
  result += d_dim * (d_dim - 1.0) * M_LN_SQRT_PI;
  for (int i = 1; i <= dim; ++i) result += 2.0 * lgammafn((degreesOfFreedom + 1.0 - (double) i) / 2.0);
  
  return(result);
}
