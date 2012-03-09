#ifndef BLME_COVARIANCE_PRIOR_H
#define BLME_COVARIANCE_PRIOR_H

#include <R.h>
#include <Rinternals.h>

double calculateCovarianceDeviance(SEXP prior, double commonScale,
                                   const double *parameters, int factorDimension);
                                   
void setCovarianceConstraints(SEXP prior, double *boxConstraints, int factorDimension);


#endif // BLME_COVARIANCE_PRIOR_H
