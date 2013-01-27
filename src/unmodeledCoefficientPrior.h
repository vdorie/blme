#ifndef BLME_UNMODELED_COEFFICIENT_PRIOR_H
#define BLME_UNMODELED_COEFFICIENT_PRIOR_H

#include <R.h>
#include <Rdefines.h>

// expects commonScale as a variance; returns only the non-constant part
double getUnmodeledCoefficientDevianceVaryingPart(SEXP prior, double commonScale,
                                                  const double *unmodeledCoefficients, int numUnmodeledCoefs);
double getUnmodeledCoefficientDevianceConstantPart(SEXP prior, int numUnmodeledCoefs);

double getUnmodeledCoefficientDensityExponentialPart(SEXP prior, const double* unmodeledCoefficients, int numUnmodeledCoefs);


// common scale has form (sigma^2)^-(df/2) * exp(-0.5 * sigma^-2 * stuff) in the objective function
// adjustment to power of common scale due to prior in terms of "degrees of freedom"
double getUnmodeledCoefficientPriorCommonScaleDegreesOfFreedom(SEXP prior, int numUnmodeledCoefs);

void addGaussianContributionToDenseBlock(SEXP regression, double *lowerRightBlock, double commonScale);

#endif // BLME_UNMODELED_COEFFICIENT_PRIOR_H
