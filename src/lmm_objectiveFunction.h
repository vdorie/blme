#ifndef BLME_LMM_OBJECTIVEFUNCTION_H
#define BLME_LMM_OBJECTIVEFUNCTION_H

#include <R.h>
#include <Rdefines.h>

#include "cache.h"

void weightMatrices(SEXP regression, MERCache* cache);
void updateMatrixFactorizations(SEXP regression, MERCache* cache);
void calculateFirstHalfOfProjections(SEXP regression, MERCache* cache);

void calculateSecondHalfOfProjections(SEXP regression, MERCache*cache);
void calculateSumsOfSquares(SEXP regression, MERCache* cache);

#endif // BLME_LMM_OBJECTIVEFUNCTION_H
