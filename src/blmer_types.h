#ifndef BLME_BLMER_TYPES_H
#define BLME_BLMER_TYPES_H

#define PRIOR_TYPE_SLOT(x)      (((priorType_t *) INTEGER(GET_SLOT((x), blme_prior_typeSym)))[0])
#define PRIOR_FAMILIES_SLOT(x) ((priorFamily_t *) INTEGER(GET_SLOT((x), blme_prior_familiesSym)))
#define PRIOR_SCALES_SLOT(x)             ((int *) INTEGER(GET_SLOT((x), blme_prior_scalesSym)))
#define PRIOR_HYPERPARAMETERS_SLOT(x) REAL(GET_SLOT((x), blme_prior_hyperparametersSym))

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

// not an enum as much as a bit mask
#define PRIOR_SCALE_POSTERIOR_MASK 0x00000001
#define PRIOR_SCALE_COMMON_MASK    0x00000002

#define getPosteriorScaleBit(x) ((priorPosteriorScale_t) ((x) & PRIOR_SCALE_POSTERIOR_MASK))
#define getCommonScaleBit(x)    ((priorCommonScale_t) (((x) & PRIOR_SCALE_COMMON_MASK) >> 1))
// for univariates situations, whether or not posterior is of sd or var
// has no effect on multivariate
typedef enum {
  PRIOR_POSTERIOR_SCALE_SD  = 0,  // must be first
  PRIOR_POSTERIOR_SCALE_VAR = 1,
  
  PRIOR_POSTERIOR_SCALE_END,
} priorPosteriorScale_t;

// whether or not the prior applies to the variance with the common scale
// multiplied in or not
typedef enum {
  PRIOR_COMMON_SCALE_FALSE = 0, // must be first
  PRIOR_COMMON_SCALE_TRUE  = 1,
  
  PRIOR_COMMON_SCALE_END,
} priorCommonScale_t;

#endif /* BLME_BLMER_H */
