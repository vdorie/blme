#ifndef BLME_BLMER_H
#define BLME_BLMER_H

#include <R.h>
#include <Rdefines.h>

#include "parameters.h"
#include "cache.h" // typedefs MERCache
#include "blmer_types.h"

int getNumParametersForParameterization(SEXP regression, parameterizationType_t type);

void setBoxConstraints(SEXP regression, double* boxConstraints);
void convertOptimizationParametersToConvergence(SEXP regression, const double* source, double* target);
void copyParametersFromRegression(SEXP regression, double* parameters);
void copyParametersIntoRegression(SEXP regression, const double* parameters);

int priorsProhibitParameters(SEXP regression, double* parameters);

double getPriorPenalty(SEXP regression, MERCache* cache, const double* parameters);

void guaranteeValidPrior(SEXP regression); // installs an empty prior if none is already in model,

int parametersIncludeUnmodeledCoefs(SEXP regression);
int parametersIncludeCommonScale(SEXP regression);
// int canProfileCommonScale(SEXP regression);
// int commonScaleRequiresOptimization(SEXP regression);

// the following exposes our enumerations to R in the form an ordered vector of strings
SEXP bmer_getTypeEnumeration();
SEXP bmer_getFamilyEnumeration();
SEXP bmer_getPosteriorScaleEnumeration();
SEXP bmer_getCommonScaleEnumeration();
SEXP bmer_getScaleInt(SEXP posteriorScale, SEXP commonScale);

SEXP bmer_getObjectiveFunction(SEXP regression);
SEXP bmer_getObjectiveFunctionForFixedCommonScale(SEXP regression);
SEXP bmer_getPriorPenalty(SEXP regression);
SEXP bmer_getCommonScaleDerivatives(SEXP regression);
SEXP bmer_getOptimalCommonScale(SEXP regression);

void printAllPriors(SEXP regression);

#endif /* BLME_BLMER_H */
