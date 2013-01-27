#ifndef BLME_COVARIANCE_PRIOR_H
#define BLME_COVARIANCE_PRIOR_H

#include <R.h>
#include <Rinternals.h>

#include "cache.h" // typedefs the MERCache
#include "blmer.h"

double getCovarianceDeviance(SEXP regression, const double* parameters);

// these two functions split up the deviance into the normalizing constant
// and whatever depends on the parameter
double getCovarianceDevianceConstantPart(SEXP regression);
double getCovarianceDevianceVaryingPart(SEXP regression, const double* parameters, int* numParametersUsed);

int getNumCovarianceParametersForPrior(priorType_t type, int dim);

// the following two are largely for lmms

// the objective func can look like
//   det(sigma * Sigma) * exp(-0.5 * stuff / sigma^2 - stuff / sigma ...)
//
// this figures out what part of "stuff" is due to the covariance prior and
// sticks them into the cache
void setCovariancePriorCommonScaleExponentialVaryingParts(SEXP regression, MERCache* cache, const double* parameters);

// get the rest, aka the determinants for all and exponential parts for priors not on the common scale
double getCovarianceDevianceCommonScaleFreeVaryingPart(SEXP regression, const double* parameters);

// common scale has form (sigma^2)^-(df/2) * exp(-0.5 * sigma^-2 * stuff) in the objective function
// adjustment to power of common scale due to prior
double getCovariancePriorCommonScaleDegreesOfFreedom(SEXP regression);
                                   
void setCovarianceConstraints(SEXP prior, double* boxConstraints, int factorDimension);

#endif // BLME_COVARIANCE_PRIOR_H
