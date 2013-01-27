#ifndef BLME_DIST_WISHART_H
#define BLME_DIST_WISHART_H

// w is the argument. it should be a dim x dim positive definite matrix,
// as should the scale parameter.
// useLog is a boolean for the scale of the result.
double dwish(const double *w, int dim, double degreesOfFreedom,
             const double *scale, int useLog);

double diwish(const double *w, int dim, double degreesOfFreedom,
              const double *inverseScale, int useLog);

// for less-R based, more internal calculations using pre-computed pieces
//
// acts like above but returns -2 * log of density and splits the result into a
// constant part and a part that depends on the random variable
// accepts input in "ST" matrix storage
double wishartDevianceVaryingPart(const double* st, int dim, double degreesOfFreedom, const double* scaleInverse);
double wishartDevianceConstantPart(int dim, double degreesOfFreedom, double logScaleDeterminant);
double wishartDevianceVaryingPartWithScale(const double* st, int dim, double degreesOfFreedom, const double* scaleInverse, double scale); // "scale" is a scalar

// further chops the varying part into an exponential and determinant term
double wishartDevianceVaryingDeterminantPart(const double* st, int dim, double degreesOfFreedom);
double wishartDevianceVaryingExponentialPart(const double* st, int dim, const double* scaleInverse);

double inverseWishartDevianceVaryingPart(const double* st, int dim, double degreesOfFreedom, const double *inverseScale);
double inverseWishartDevianceConstantPart(int dim, double degreesOfFreedom, double logInverseScaleDeterminant);
double inverseWishartDevianceVaryingPartWithScale(const double* st, int dim, double degreesOfFreedom, const double* scaleInverse, double scale);

double inverseWishartDevianceVaryingDeterminantPart(const double* st, int dim, double degreesOfFreedom);
double inverseWishartDevianceVaryingExponentialPart(const double* st, int dim, const double* scaleInverse);

#endif /* BLME_DIST_WISHART_H */
