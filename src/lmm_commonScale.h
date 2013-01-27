#ifndef BLME_LMM_COMMONSCALE_H
#define BLME_LMM_COMMONSCALE_H

#include <R.h>
#include <Rinternals.h>

#include "cache.h"

#define COMMON_SCALE_OPTIMIZATION_TOLERANCE 1.0e-10

// All of the below require the cache to filled out in some way (not just constructed).
// More or less, they only make sense in the context of lmmGetObjectiveFunction in
// objectiveFunction.c.

// optimize will do so numerically if necessary, otherwise will profile
void optimizeCommonScale(SEXP regression, MERCache* cache);
// profile is just profile; if called when inappropriate, does nothing
void profileCommonScale(SEXP regression, MERCache* cache);

// below are mostly for testing purposes and theory work
double performOneStepOfNewtonsMethodForCommonScale(SEXP regression, MERCache* cache);
void getCommonScaleDerivatives(SEXP regression, MERCache* cache, double* firstDerivative, double* secondDerivative);

#endif // BLME_LMM_COMMONSCALE_H
