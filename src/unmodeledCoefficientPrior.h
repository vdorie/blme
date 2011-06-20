#ifndef BLME_UNMODELED_COEFFICIENT_PRIOR_H
#define BLME_UNMODELED_COEFFICIENT_PRIOR_H

#include <R.h>
#include <Rdefines.h>

double
calculateUnmodeledCoefficientDeviance(SEXP prior, double commonScale,
                                      const double *unmodeledCoefficients, int numUnmodeledCoefs);

void addGaussianContributionToDenseBlock(SEXP regression, double *lowerRightBlock);

#endif // BLME_UNMODELED_COEFFICIENT_PRIOR_H
