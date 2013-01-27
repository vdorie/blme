#ifndef BLME_LMM_PRIORS_H
#define BLME_LMM_PRIORS_H

#include <R.h>
#include <Rinternals.h>

#include "cache.h"

void setPriorConstantContributionsToCommonScale(SEXP regression, MERCache* cache);
void setPriorVaryingContributionsToCommonScale(SEXP regression, MERCache* cache, const double* parameters);

double getPriorDevianceCommonScaleFreeVaryingPart(SEXP regression, const double* parameters);
double getPriorDevianceConstantPart(SEXP regression);

// not for typical use, calculates everything form scratch and skips over most optimizations
double getPriorDeviance(SEXP regression, MERCache* cache, const double* parameters);

#endif // BLME_LMM_PRIORS_H
