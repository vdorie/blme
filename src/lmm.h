#ifndef BLME_LMM_H
#define BLME_LMM_H

#include <R.h>
#include <Rinternals.h>

#include "cache.h"

double lmmGetObjectiveFunction(SEXP regression, MERCache* cache, const double* parameters);
// when calculating the objective function, a lot of stuff gets put into the cache that
// might be of interest. The below cleans that stuff up and sticks it into the MER object
void lmmFinalizeOptimization(SEXP regression, MERCache* cache);

double lmmGetObjectiveFunctionForFixedCommonScale(SEXP regression, MERCache* cache, const double* parameters);
double lmmGetPriorPenalty(SEXP regression, MERCache* cache, const double* parameters);

// legacy lmer-ey thing
void lmmCalculateProjectionsForSingleArgumentAnova(SEXP regression, double *modeledCoefProjection, double *unmodeledCoefProjection);


// X'X - Rzx'Rzx. externally available for sim(). also used in objective func
void computeDowndatedDenseCrossproduct(SEXP regression, MERCache* cache, double* target);

MERCache* createLMMCache(SEXP regression);
void deleteLMMCache(MERCache* cache);


// for testing purposes and theoretical work; probably best if not used, hence the '_' prefix
void _getCommonScaleDerivatives(SEXP regression, MERCache* cache, const double* parameters, double* firstDerivative, double* secondDerivative);
double _getOptimalCommonScale(SEXP regression, MERCache* cache, const double* parameters);

#endif // BLME_LMM_H
