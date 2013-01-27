#include <Rdefines.h>

#include <Rmath.h>               /* gamma functions and finite/NaN check */
#include <R_ext/Lapack.h>        /* for Lapack and BLAS */

#include "dist_gamma.h"

double gammaDevianceConstantPart(double shape, double rate)
{
  if (shape == 0.0 || rate == 0.0) return(0.0);
  if (shape <  0.0 || rate <  0.0) return(R_NaN);
  
  return(-2.0 * (shape * log(rate) - lgammafn(shape)));
}

double gammaDevianceVaryingPart(double x, double shape, double rate)
{
  if (x == 0.0) return(shape > 1.0 ? R_PosInf : R_NegInf);
  if (x <  0.0) return(R_NaN);
  
  return(-2.0 * ((shape - 1.0) * log(x) - rate * x));
}

double gammaDevianceVaryingPolynomialPart(double x, double shape)
{
  if (x == 0.0) return(shape > 1.0 ? R_PosInf : R_NegInf);
  if (x <  0.0) return(R_NaN);
  
  return(-2.0 * (shape - 1.0) * log(x));
}

double gammaDevianceVaryingExponentialPart(double x, double rate)
{
  if (x < 0.0) return(R_NaN);
  
  return( 2.0 * rate * x);
}


double inverseGammaDevianceConstantPart(double shape, double scale)
{
  if (shape == 0.0 || scale == 0.0) return(0.0);
  if (shape <  0.0 || scale <  0.0) return(R_NaN);
  
  return(-2.0 * (shape * log(scale) - lgammafn(shape)));
}

double inverseGammaDevianceVaryingPart(double x, double shape, double scale)
{
  if (x == 0.0) return(scale > 0.0 ? R_PosInf : R_NegInf);
  if (x <  0.0) return(R_NaN);
  
  return(2.0 * ((shape + 1.0) * log(x) + scale / x));
}

double inverseGammaDevianceVaryingPolynomialPart(double x, double shape)
{
  if (x == 0.0) return(R_NegInf);
  if (x <  0.0) return(R_NaN);
  
  return(2.0 * (shape + 1.0) * log(x));
}

double inverseGammaDevianceVaryingExponentialPart(double x, double scale)
{ 
  if (x == 0.0) {
    if (scale == 0.0) return(R_NaN); // 0 / 0, more or less
    return(R_PosInf); // limit of -2 x / 0
  }
  if (x <  0.0) return(R_NaN);
  
  return(2.0 * scale / x);
}
