#ifndef BLME_LMM_H
#define BLME_LMM_H

#include <R.h>
#include <Rinternals.h>

typedef struct _MERCache MERCache; // defined in lmm.c

MERCache *createLMMCache(SEXP regression);
void deleteLMMCache(MERCache *cache);

double lmmCalculateDeviance(SEXP regression, MERCache *cache);
double lmmApproximateDeviance(SEXP regression, MERCache *cache);

void updateWeights(SEXP regression, MERCache *cache);
void updateAugmentedDesignMatrixFactorizations(SEXP regression, MERCache *cache);

// We project the responses *into* the normalized column space of the modeled
// and unmodeled design matrices column space.
//
// If (X'X)^-1 X' y is the usual projection into, by replacing X with X* = X L^-T,
// LL' = X'X, we "normalize" the column space, i.e.:
//   X*'X* = L^-1 X'X L^-T = L^-1 LL' L^-T = I
//
// In general, we are not interested in the normalized projections except for
// the fact that they have same length as the projections *onto* the column
// space of the design matrices. See the accompanying pdf for the details of
// optimization in the linear case.
void 
calculateProjections(SEXP regression, MERCache *cache);
void
calculatePenalizedWeightedResidualSumOfSquaresFromProjections(SEXP regression, MERCache *cache);
// we rotate back out of the normalized space to obtain the projection
// onto the regular column space
void
rotateProjections(SEXP regression, MERCache *cache);

// X'X - Rzx'Rzx'. externally available for sim()
void computeDowndatedDenseCrossproduct(SEXP regression, MERCache* cache, double* target);

void profileCommonScale(SEXP regression, MERCache* cache);
void updateDeviance(SEXP regression, MERCache* cache);

// returns the new common scale
double performOneStepOfNewtonsMethodForCommonScale(SEXP regression, MERCache *cache);


void calculateProjectionsForSingleArgumentAnova(SEXP regression,
                                                double *modeledCoefProjection,
                                                double *unmodeledCoefProjection);

                                                
#endif // BLME_LMM_H
