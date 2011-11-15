#ifndef BLME_BLMER_H
#define BLME_BLMER_H

#include <R.h>
#include <Rdefines.h>
#include "parameters.h"

#define PRIOR_TYPE_SLOT(x) (((priorType_t *) INTEGER(GET_SLOT((x), blme_prior_typeSym)))[0])
#define PRIOR_FAMILIES_SLOT(x) ((priorFamily_t *) INTEGER(GET_SLOT((x), blme_prior_familiesSym)))
#define PRIOR_SCALES_SLOT(x) ((priorScale_t *) INTEGER(GET_SLOT((x), blme_prior_scalesSym)))
#define PRIOR_HYPERPARAMETERS_SLOT(x) REAL(GET_SLOT((x), blme_prior_hyperparametersSym))

int getNumParametersForParameterization(SEXP regression, parameterizationType_t type);

void convertOptimizationParametersToConvergence(SEXP regression, const double *source, double *target);
void initializeOptimizationParameters(SEXP regression, double *parameters);
void updateRegressionWithParameters(SEXP regression, const double *parameters);
void setBoxConstraints(SEXP regression, double *boxConstraints);

int isAtBoundary(SEXP regression, double *parameters);

double calculatePriorPenalty(SEXP regression, double *parameters);

void guaranteeValidPrior(SEXP regression); // installs an empty prior if none is already in model,

int parametersIncludeUnmodeledCoefs(SEXP regression);
int parametersIncludeCommonScale(SEXP regression);
int canProfileCommonScale(SEXP regression);
int commonScaleRequiresOptimization(SEXP regression);

// the following exposes our enumerations to R in the form an ordered vector of strings
SEXP bmer_getTypeEnumeration();
SEXP bmer_getFamilyEnumeration();
SEXP bmer_getScaleEnumeration();

SEXP bmer_calculatePriorPenalty(SEXP regression);

void printAllPriors(SEXP regression);

// enumerations
typedef enum {
  PRIOR_TYPE_NONE = 0,    // must be first, p(x) \propto 1, i.e. just return the current estimate
  PRIOR_TYPE_DIRECT,      // no decomposition, p(Sigma) = whatever is named
  PRIOR_TYPE_CORRELATION, // work on a SRS = sd*correlation*sd decomposition
  PRIOR_TYPE_SPECTRAL,    // perform a QLambdaQ' = orthogonal*diagonal(eigenvalues)*orthogonal' decomposition
  
  PRIOR_TYPE_END
} priorType_t;

typedef enum {
  PRIOR_FAMILY_FLAT = 0,     // flat must be first
  
  PRIOR_FAMILY_GAMMA,        // priors on var/cov
  PRIOR_FAMILY_INVGAMMA,
  PRIOR_FAMILY_WISHART,
  PRIOR_FAMILY_INVWISHART,
  
  PRIOR_FAMILY_GAUSSIAN,     // priors on fixef
  PRIOR_FAMILY_MVT,
  
  PRIOR_FAMILY_POINT,        // priors on common scale
  
  PRIOR_FAMILY_END
} priorFamily_t;

typedef enum {
  PRIOR_SCALE_SD = 0,     // sd must be first; relates to variances
  PRIOR_SCALE_VARIANCE,
  
  PRIOR_SCALE_ABSOLUTE,   // relates to variance of unmodeled coef prior
  PRIOR_SCALE_COMMON,
  
  PRIOR_SCALE_END
} priorScale_t;

#endif /* BLME_BLMER_H */
