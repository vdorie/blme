#ifndef LME4_LMER_H
#define LME4_LMER_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>

// #include "Syms.h"

/** positions in the deviance vector */
enum devP {
  ML_POS=0,			/**<Maximum likelihood estimation criterion  */
  REML_POS,			/**<REML criterion */
  ldL2_POS,			/**<2*log-determinant of L */
  ldRX2_POS,			/**<2*log-determinant of RX */
  sigmaML_POS,		/**<current ML estimate of sigma */
  sigmaREML_POS,		/**<current REML estimate of sigma */
  pwrss_POS,			/**<penalized weighted residual sum of squares */
  disc_POS,			/**<discrepancy */
  usqr_POS,			/**<squared length of u */
  wrss_POS,			/**<weighted residual sum of squares  */
  dev_POS,			/**<deviance - defined for quasi families  */
  llik_POS,			/**<log-likelihood - undefined for quasi families  */
  NULLdev_POS			/**<null deviance */
};
/** positions in the dims vector */
enum dimP {
  nt_POS=0,			/**<number of terms in random effects */
  n_POS,			/**<number of observations */
  p_POS,			/**<number of fixed-effects parameters */
  q_POS,			/**<number of random effects */
  s_POS,			/**<number of variables in h (1 unless nonlinear) */
  np_POS,			/**<total number of parameters for T and S */
  LMM_POS,			/**<is the model a linear mixed model? */
  isREML_POS,			/**<indicator of REML estimation */
  fTyp_POS,			/**<family type for generalized model */
  lTyp_POS,			/**<link type for generalized model */
  vTyp_POS,			/**<variance type for generalized model */
  nest_POS,			/**<indicator of nested grouping factors */
  useSc_POS,			/**<does the family use a separate scale parameter */
  nAGQ_POS,			/**<number of adaptive Gauss-Hermite quadrature pts */
  verb_POS,			/**<verbose output in mer_optimize? */
  mxit_POS,			/**<maximum # of iterations in mer_optimize */
  mxfn_POS,			/**<maximum # of function evaluations in mer_optimize */
  cvg_POS,			/**<convergence indictor from port optimization  */
};

SEXP lme4_ghq(SEXP np);
SEXP mer_ST_getPars(SEXP x);
SEXP mer_ST_initialize(SEXP ST, SEXP Gp, SEXP Zt);
SEXP mer_ST_setPars(SEXP x, SEXP pars);
SEXP mer_update_dev(SEXP x);
SEXP mer_create_L(SEXP CmP);
SEXP mer_optimize(SEXP x);
SEXP mer_update_ranef(SEXP x);
SEXP mer_update_mu(SEXP x);
SEXP mer_validate(SEXP x);

#endif /* LME4_LMER_H */
