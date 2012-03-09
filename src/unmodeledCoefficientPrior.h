#ifndef BLME_UNMODELED_COEFFICIENT_PRIOR_H
#define BLME_UNMODELED_COEFFICIENT_PRIOR_H

#include <R.h>
#include <Rdefines.h>

#include "lmm.h" // for MERCache

double
calculateUnmodeledCoefficientDeviance(SEXP prior, double commonScale,
                                      const double *unmodeledCoefficients, int numUnmodeledCoefs);

void addGaussianContributionToDenseBlock(SEXP regression, double *lowerRightBlock, double commonScale);


// derivatives of objective function w/r/t common scale
void getDerivatives(SEXP regression, MERCache* cache,
                    double* firstDerivative, double* secondDerivative);

// externally callable version of above
SEXP bmer_getDerivatives(SEXP regression);

#endif // BLME_UNMODELED_COEFFICIENT_PRIOR_H
