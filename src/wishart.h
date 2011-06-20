#ifndef BLME_WISHART_H
#define BLME_WISHART_H

// w is the argument. it should be a dim x dim positive definite matrix,
// as should the scale parameter.
// useLog is a boolean for the scale of result.
double dwish(const double *w, int dim, double degreesOfFreedom,
             const double *scale, int useLog);

// same as above, but the user can provide a few intermediate
// steps in the calculation
//
// note, this uses the scale matrix inverse
// conversely, the inverse wishart is parameterized by the inverse of a scale matrix
// the distinction is that in the first, the user supplies the inverse of a parameter,
// while in the second, it is the parameter itself that is supplied
double dwish2(const double *w, int dim, double degreesOfFreedom,
             double logScaleDeterminant, const double *scaleInverse, int useLog);

double diwish(const double *w, int dim, double degreesOfFreedom,
              const double *inverseScale, int useLog);
              
double diwish2(const double *wInverse, int dim, double degreesOfFreedom,
               double logInverseScaleDeterminant, const double *inverseScale, int useLog);

#endif /* BLME_WISHART_H */
