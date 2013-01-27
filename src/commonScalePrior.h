#ifndef BLME_COMMON_SCALE_PRIOR_H
#define BLME_COMMON_SCALE_PRIOR_H

#include <R.h>
#include <Rdefines.h>

#include "cache.h" // typedefs MERCache

// expects commonScale as a variance; returns only the non-constant part
double getCommonScaleDevianceVaryingPart(SEXP prior, double commonScale);
// gets the constant part. Add this to the above to get the full density
double getCommonScaleDevianceConstantPart(SEXP prior);

// common scale has form (sigma^2)^-(df/2) * exp(-0.5 * sigma^-2 * stuff) in the objective function
// adjustment to power of common scale due to prior
double getCommonScalePriorDegreesOfFreedom(SEXP prior);
// adjustment to exponent
void setCommonScalePriorExponentialConstantPart(SEXP prior, MERCache* cache);


#endif // BLME_COMMON_SCALE_PRIOR_H
