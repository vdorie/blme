/* GLOSSARY OF TERMS/ABBREVIATIONS
 * -------------------------------
 *
 * Abbreviations:
 *   coef(s)   - coefficient(s)
 **/

/* Further notes on notation
 * -------------------------
 *
 * The "augmented" design matrix arises when writing the joint distribution
 * of the modeled coefficients and the observations in a linear, multilevel
 * model. Specifially, the exponential term can be written as the norm of:
 *
 *   [ ZLambda X ] [ u    ] - [ y ]
 *   [ I       0 ] [ beta ] - [ 0 ]
 *
 * This takes the form of a simple linear regression, and the "design matrix"
 * of the vector of modeled and unmodeled coefficients we call augmented.
 *
 * We compute the cholesky factorization of the crossproduct of this matrix.
 * In terms of the original lmer code, we have that
 *
 *   [ L      0 ] [ L' RZX ] = [ Lambda'Z'ZLambda + I Lambda'Z'X ]
 *   [ RZX' RX' ] [ 0   RX ] = [ X'ZLambda            X'X        ]
 *
 * Returning to the original design matrix, we have that the first block-column
 * [ Lambda'Z' I ]' is sparse, while the second [ X' 0 ]' is dense.
 * Consequently, the upper left block of the crossproduct is again sparse, while
 * the other factors are dense. In the code "L", "RZX", and "RX" thus called the
 * sparse left factor, the off-diagonal right factor, and the dense right factor.
 */
 
#include <Rmath.h>		 /* for dnorm5, etc. */
#include <R_ext/Lapack.h>        /* for Lapack (dpotrf, etc.) and BLAS */
#include <R_ext/stats_package.h> /* for S_nlminb_iterate */
#include <Matrix.h>		 /* for cholmod functions */

#include "lmer.h"
#include "Syms.h"
#include "blmer.h"
#include "util.h"
#include "lmm.h"
#include "lmer_common.h"

#include "common_inlines.h"

extern	       // cholmod_common struct initialized in R_init_lme4 
cholmod_common cholmodCommon;

#ifdef ENABLE_NLS		// Allow for translation of error messages
#include <libintl.h>
#define _(String) dgettext ("lme4", String)
#else
#define _(String) (String)
#endif

/* Constants */

#ifndef BUF_SIZE
// size of buffer for an error message
#define BUF_SIZE 127
#endif	

// Maximum number of iterations in update_u
#define CM_MAXITER  300
// Tolerance level for convergence criterion in update_u
#define CM_TOL      1e-10
// Minimum step factor in update_u
#define CM_SMIN     1e-5

// precision and maximum number of iterations used in GHQ
#define GHQ_EPS    1e-15
#define GHQ_MAXIT  40

#define LTHRESH     30.
#define MLTHRESH   -30.

static double MPTHRESH = 0;
static double PTHRESH = 0;
static const double INVEPS = 1/DOUBLE_EPS;

/* In-line functions */


/**
 * Evaluate y * log(y/mu) with the correct limiting value at y = 0.
 *
 * @param y 
 * @param mu
 *
 * @return y * log(y/mu) for y > 0, 0 for y == 0.
 */
static R_INLINE double y_log_y(double y, double mu)
{
  return (y) ? (y * log(y/mu)) : 0;
}

/* Stand-alone utility functions (sorted by function name) */

/**
 * Check that slot sym of object x is a numeric matrix of dimension nr
 * by nc.
 *
 * @param buf character buffer of length nb + 1
 * @param nb number of writable positions in the buffer
 * @param x pointer to an mer object
 * @param sym name (symbol, actually) of the slot to check
 * @param nr expected number of rows
 * @param nc expected number of columns
 *
 * @return 0 for success, number of characters written to buf on failure
 */
static int chkDims(char *buf, int nb, SEXP x, SEXP sym, int nr, int nc)
{
  SEXP MP = GET_SLOT(x, sym);
  int *dm = isMatrix(MP) ? INTEGER(getAttrib(MP, R_DimSymbol)) : (int*) NULL;
  if (!dm || !isReal(MP) || dm[0] != nr || dm[1] != nc)
    return snprintf(buf, BUF_SIZE,
                    _("Slot %s must be a numeric matrix of size %d by %d"),
                    CHAR(PRINTNAME(sym)), nr, nc);
  return 0;
}

/** Check that the length of the sym slot in x is len or, possibly, zero.
 *
 * @param buf character buffer of length nb + 1
 * @param nb number of writable positions in the buffer
 * @param x pointer to an mer object
 * @param sym name (symbol, actually) of the slot to check
 * @param len expected length
 * @param zerok is a length of zero allowed?
 *
 * @return 0 for success, number of characters written to buf on failure
 
 */
static int chkLen(char *buf, int nb, SEXP x, SEXP sym, int len, int zerok)
{
  int ll;
  
  if (!(ll = LENGTH(GET_SLOT(x, sym))) == len && zerok && !ll)
    return snprintf(buf, BUF_SIZE, _("Slot %s must have length %d."),
                    CHAR(PRINTNAME(sym)), len);
  return 0;
}

/**
 * Evaluate the sum of the deviance residuals for a GLM family
 *
 * @param mu pointer to the mu vector
 * @param pWt pointer to the vector of prior weights (NULL for
 *            constant weights)
 * @param y pointer to the response vector
 * @param n length of mu and y
 * @param vTyp type of variance function: the 1-based index into
 *        c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")
 * @param ans pointer to vector of partial sums
 * @param fac indices associating observations with partial sums 
 *            (may be (int*)NULL)
 * @return the sum of the deviance residuals
 */
static double*
lme4_devResid(const double* mu, const double* pWt, const double* y,
              int n, int vTyp, double *ans, int *fac)
{
  for (int i = 0; i < n; i++) {
    double mu_i = mu[i];
    double weight_i = pWt ? pWt[i] : 1.0;
    double y_i = y[i];
    double residual_i = y_i - mu_i;
    
    int a_i = fac ? (fac[i] - 1) : 0;
    switch(vTyp) {
      case 1:			/* constant variance */
        ans[a_i] += weight_i * residual_i * residual_i;
        break;
      case 2:			/* mu(1-mu) variance */
        ans[a_i] += 2 * weight_i *
          (y_log_y(y_i, mu_i) + y_log_y(1 - y_i, 1 - mu_i));
        break;
      case 3:			/* mu variance */
        ans[a_i] += 2 * weight_i * (y_log_y(y_i, mu_i) - (y_i - mu_i));
        break;
      case 4:			/* mu^2 variance */
        ans[a_i] += 2 * weight_i * (y_log_y(y_i, mu_i) - (y_i - mu_i)/mu_i);
        break;
      case 5:			/* mu^3 variance */
        ans[a_i] += weight_i * (residual_i * residual_i)/(y_i * mu_i * mu_i);
        break;
      default:
        error(_("Unknown vTyp value %d"), vTyp);
    }
  }
  return ans;
}

/**
 * Evaluate the derivative d mu/d eta for the GLM link function of
 * type lTyp
 *
 * @param mu pointer to the mu vector
 * @param muEta pointer to the muEta vector
 * @param eta pointer to the eta vector
 * @param n length of mu, muEta and eta
 * @param lTyp type of link: the 1-based index into
 *        c("logit", "probit", "cauchit", "cloglog", "identity",
 *          "log", "sqrt", "1/mu^2", "inverse")
 */
static void
lme4_muEta(double* mu, double* muEta, const double* eta, int n, int lTyp)
{
  for (int i = 0; i < n; i++) { /* apply the generalized linear part */
    double etai = eta[i], tmp, t2;
    switch(lTyp) {
      case 1:		/* logit */
        tmp = (etai < MLTHRESH) ? DOUBLE_EPS : ((etai > LTHRESH) ?
                                                INVEPS : exp(etai));
        mu[i] = tmp/(1 + tmp);
        muEta[i] = mu[i] * (1 - mu[i]);
        break;
      case 2:		/* probit */
        if (!MPTHRESH) {
          MPTHRESH = qnorm5(DOUBLE_EPS, 0, 1, 1, 0);
          PTHRESH = -MPTHRESH;
        }
        mu[i] = (etai < MPTHRESH) ? DOUBLE_EPS :
          ((etai > PTHRESH) ? 1 - DOUBLE_EPS :
           pnorm5(etai, 0, 1, 1, 0));
        tmp = dnorm4(eta[i], 0, 1, 0);
        muEta[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
        break;
      case 3:		/* cauchit */
        error(_("cauchit link not yet coded"));
        break;
      case 4:		/* cloglog */
        tmp = (etai < MLTHRESH) ? DOUBLE_EPS : ((etai > LTHRESH) ?
                                                INVEPS : exp(etai));
        t2 = -expm1(-tmp);
        mu[i] = (t2 < DOUBLE_EPS) ? DOUBLE_EPS : t2;
        muEta[i] = tmp * exp(-tmp);
        break;
      case 5:		/* identity */
        mu[i] = etai;
        muEta[i] = 1.;
        break;
      case 6:		/* log */
        tmp = exp(etai);
        muEta[i] = mu[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
        break;
      case 7:		/* sqrt */
        mu[i] = etai * etai;
        muEta[i] = 2 * etai;
        break;
      case 8:		/* 1/mu^2 */
        mu[i] = sqrt(etai);
        muEta[i] = 1/(2*sqrt(etai));
        break;
      case 9:		/* inverse */
        mu[i] = 1/etai;
        muEta[i] = -1/(etai*etai);
        break;
      default:
        error(_("General form of glmer_linkinv not yet written"));
    }
  }
}

/**
 * Evaluate the GLM variance function of type vTyp given mu
 *
 * @param var pointer to the variance vector
 * @param mu pointer to the mu vector
 * @param n length of var and mu
 * @param vTyp type of variance function: the 1-based index into
 *        c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")
 *
 */
static void
lme4_varFunc(double* var, const double* mu, int n, int vTyp)
{
  for (int i = 0; i < n; i++) {
    double mui = mu[i];
    switch(vTyp) {
      case 1:			/* constant variance */
        var[i] = 1.;
        break;
      case 2:			/* mu(1-mu) variance */
        if (mui <= 0 || mui >= 1)
          error(_("mu[i] must be in the range (0,1): mu = %g, i = %d"),
                mu, i);
        var[i] = mui * (1 - mui);
        break;
      case 3:			/* mu variance */
        if (mui <= 0)
          error(_("mu[i] must be positive: mu = %g, i = %d"), mu, i);
        var[i] = mui;
        break;
      case 4:			/* mu^2 variance */
        if (mui <= 0)
          error(_("mu[i] must be positive: mu = %g, i = %d"), mu, i);
        var[i] = mui * mui;
        break;
      case 5:			/* mu^3 variance */
        if (mui <= 0)
          error(_("mu[i] must be positive: mu = %g, i = %d"), mu, i);
        var[i] = mui * mui * mui;
        break;
      default:
        error(_("Unknown vTyp value %d"), vTyp);
    }
  }
}

/**
 * Extract the parameters from ST list
 *
 * @param x an mer object
 * @param pars vector of the appropriate length
 *
 * @return pointer to the parameter vector
 */
static double *ST_getPars(SEXP x, double *pars)
{
  SEXP ST = GET_SLOT(x, lme4_STSym);
  int nt = LENGTH(ST), pos = 0;
  for (int i = 0; i < nt; i++) {
    SEXP STi = VECTOR_ELT(ST, i);
    double *st = REAL(STi);
    int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
    int ncp1 = nci + 1;
    
    for (int j = 0; j < nci; j++)
	    pars[pos++] = st[j * ncp1];
    for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
        pars[pos++] = st[k + j * nci];
  }
  return pars;
}

/**
 * Populate the st, numColumnsPerFactor and numLevelsPerFactor arrays.  Return the maximum element of nc.
 *
 * @param ST pointer to a list (length nt) of matrices
 * @param Gp group pointers (length nt + 1)
 * @param st length nt array of (double*) pointers to be filled with
 * pointers to the contents of the matrices in ST.  Not used if NULL.
 * @param nc length nt array to be filled with the number of columns
 * @param nlev length nt array to be filled with the number of
 *        levels of the grouping factor for each term
 * 
 * @return maximum element of nc
 */
static int
extractSTParametersAndStructure(const SEXP stExpression, const int *sparseRowForFactor,
                                double **stMatrices, int *factorDimensions,
                                int *numGroupsPerFactor)
{
  SparseMatrixStructure matrixStructure;
  matrixStructure.factorDimensions = factorDimensions;
  matrixStructure.numGroupsPerFactor = numGroupsPerFactor;
  
  getSparseContentAndStructure(stExpression, sparseRowForFactor, stMatrices, &matrixStructure);
  
  return(matrixStructure.maxFactorDimension);
}

/**
 * Determine the index in theta_S corresponding to index i in the u vector
 *
 * @param i index in u vector (0 <= i < q)
 * @param nt number of random effects terms in the model
 * @param Gp vector of group pointers into u (length nt + 1, Gp[0] == 0)
 * @param nlev vector of length nt giving the number of levels per term
 * @param spt vector of group pointers into theta_S (length nt + 1,
 *            spt[0] == 0)
 */
static R_INLINE int
theta_S_ind(int i, int nt, int *Gp, int *nlev, int *spt)
{
	int trm = getFactorForSparseRow(i, nt, Gp);
	return (spt[trm] + (i - Gp[trm]) / nlev[trm]);
}

/* Level-1 utilties that call at least one of the stand-alone utilities */

/**
 * Update the L, sqrtrWt and resid slots.  It is assumed that
 * update_mu has already been called at the current values of u and
 * the model parameters and that A has been updated.
 *
 * @param x pointer to an mer object
 *
 * @return penalized weighted residual sum of squares
 *
 */
static double update_L(SEXP x)
{
  int *dims = DIMS_SLOT(x);
  int numObservations = dims[n_POS];
  int numFixedEffects = dims[p_POS];
  int s = dims[s_POS];
  
  double *gradient = V_SLOT(x), *cx = Cx_SLOT(x);
  double *deviances = DEV_SLOT(x);
	double *residuals = RESID_SLOT(x), *mu = MU_SLOT(x), *muEta = MUETA_SLOT(x);
	double *priorWeights = PWT_SLOT(x);
  double *sqrtResidualWeight = SRWT_SLOT(x); // this is W^1/2 (in vector form), sqrt(priorWeight / var)
  double *sqrtModelRowWeight = SXWT_SLOT(x); // W^1/2 * diag(d mu / d eta)
	double *variances =  VAR_SLOT(x), *y = Y_SLOT(x), one[] = { 1.0, 0.0 };
  CHM_SP A = A_SLOT(x);
  CHM_FR L = L_SLOT(x);
  R_CheckStack();
  
  /* Update the square root weights and residuals. Re-evaluate the weighted residual sum of squares. */
  if (variances != NULL || priorWeights != NULL) {
    deviances[wrss_POS] = 0.0;
    for (int j = 0; j < numObservations; ++j) {
	    sqrtResidualWeight[j] = sqrt((priorWeights ? priorWeights[j] : 1.0) / (variances ? variances[j] : 1.0));
      
	    residuals[j] = sqrtResidualWeight[j] * (y[j] - mu[j]);
	    deviances[wrss_POS] += residuals[j] * residuals[j];
    }
  }
  if (sqrtModelRowWeight != NULL) {			/* Update sXwt - matrix of model row weights - and C */
    int *ai = (int*)A->i, *ap = (int*)A->p;
    double *ax = (double*)(A->x);
    CHM_SP C = A;
    
    for (int j = 0; j < s; ++j) { /* s == 1 unless NLMM */
	    for (int i = 0; i < numObservations; ++i) {
        int weightIndex = i + j * numObservations;
        sqrtModelRowWeight[weightIndex] =
          (sqrtResidualWeight ? sqrtResidualWeight[i] : 1.0) *
		      (muEta ? muEta[i] : 1.0) *
          (gradient ? gradient[weightIndex] : 1.0); /* gradient is NULL unless NLMM */
	    }
    }
    /* C is a scaled version of A */
    if (cx != NULL) { /* cx exists only when s = 1, which should be for all GLMM and LMM */
	    for (int j = 0; j < numObservations; j++) {
        for (int p = ap[j]; p < ap[j + 1]; p++) {
          cx[p] = ax[p] * sqrtModelRowWeight[j];
        }
      }
	    A->x = (void*)cx;
    } else {
	    int *ci, *cp;
	    C = Cm_SLOT(x);
	    R_CheckStack();
      
	    ci = (int*)C->i; cp = (int*)C->p; cx = (double*)C->x;
	    AZERO(cx, cp[numObservations]);
	    for (int j = 0; j < s; j++)
        for (int i = 0; i < numObservations; i++) {
          int ja = i + j * numObservations, pc;
          for (int pa = ap[ja]; pa < ap[ja + 1]; pa++) {
            for (pc = cp[i]; pc < cp[i + 1]; pc++) {
              if (ci[pc] == ai[pa]) break;
            }
            if (pc >= cp[i + 1])
              error(_("Structure of Cm and A are not consistent"));
            cx[pc] += ax[pa] * sqrtModelRowWeight[ja];
          }
        }
	    A = C;
    }
  }
  /* factorize AA', or CC' if C exists, where C = AW^(1/2)diag(d mu / d eta), and A = ST'Z' */
  if (!M_cholmod_factorize_p(A, one, NULL, 0 /*fsize*/, L, &cholmodCommon))
    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
          cholmodCommon.status, L->minor, L->n);

  
  deviances[ldL2_POS] = M_chm_factor_ldetL2(L);
  deviances[pwrss_POS] = deviances[usqr_POS] + deviances[wrss_POS];
  deviances[sigmaML_POS] =
    sqrt(deviances[pwrss_POS]/
         (sqrtResidualWeight ? getSumOfSquares(sqrtResidualWeight, numObservations) : (double) numObservations));
  if (gradient != NULL || muEta != NULL) { /* NLMM or GLMM */
    deviances[sigmaREML_POS] = NA_REAL;
  } else {
    deviances[sigmaREML_POS] = 
      deviances[sigmaML_POS] *
      sqrt(((double) numObservations) / ((double) (numObservations - numFixedEffects)));
  }
  
  return deviances[pwrss_POS];
}

/**
 * Update the eta, v, mu, resid and var slots according to the current
 * values of the parameters and u.  Also evaluate d[wrss_POS] using
 * the current sqrtrWt slot.  The sqrtrWt slot is changed in update_L.
 *
 * @param x pointer to an mer object
 *
 * @return penalized, weighted residual sum of squares
 */
static double update_mu(SEXP x)
{
  int *dims = DIMS_SLOT(x);
  int i1 = 1, n = dims[n_POS], p = dims[p_POS], s = dims[s_POS];
  int ns = n * s;
  double *gradient = V_SLOT(x), *d = DEV_SLOT(x), *eta = ETA_SLOT(x),
	*etaold = (double*) NULL, *mu = MU_SLOT(x),
	*muEta = MUETA_SLOT(x), *offset = OFFSET_SLOT(x),
	*srwt = SRWT_SLOT(x), *res = RESID_SLOT(x),
	*var = VAR_SLOT(x), *y = Y_SLOT(x), one[] = {1,0};
  CHM_FR L = L_SLOT(x);
  CHM_SP A = A_SLOT(x);
  CHM_DN Ptu, ceta, cu = AS_CHM_DN(GET_SLOT(x, lme4_uSym));
  R_CheckStack();
  
  if (gradient) {
    etaold = eta;
    eta = Calloc(ns, double);
  }
  /* eta := offset or eta := 0 */
  for (int i = 0; i < ns; i++) eta[i] = offset ? offset[i] : 0;
  /* eta := eta + X \beta */
  F77_CALL(dgemv)("N", &ns, &p, one, X_SLOT(x), &ns,
                  FIXEF_SLOT(x), &i1, one, eta, &i1);
  /* eta := eta + C' P' u */
  Ptu = M_cholmod_solve(CHOLMOD_Pt, L, cu, &cholmodCommon);
  ceta = N_AS_CHM_DN(eta, ns, 1);
  R_CheckStack();
  if (!M_cholmod_sdmult(A, 1 /* trans */, one, one, Ptu, ceta, &cholmodCommon))
    error(_("cholmod_sdmult error returned"));
  M_cholmod_free_dense(&Ptu, &cholmodCommon);
  
  if (gradient) {		/* evaluate the nonlinear model */
    SEXP pnames = VECTOR_ELT(GET_DIMNAMES(GET_SLOT(x, lme4_VSym)), 1),
    gg, rho = GET_SLOT(x, lme4_envSym), vv;
    int *gdims;
    
    if (!isString(pnames) || LENGTH(pnames) != s)
	    error(_("Slot V must be a matrix with %d named columns"), s);
    for (int i = 0; i < s; i++) { /* par. vals. into env. */
	    vv = findVarInFrame(rho,
                          install(CHAR(STRING_ELT(pnames, i))));
	    if (!isReal(vv) || LENGTH(vv) != n)
        error(_("Parameter %s in the environment must be a length %d numeric vector"),
              CHAR(STRING_ELT(pnames, i)), n);
	    Memcpy(REAL(vv), eta + i * n, n);
    }
    vv = PROTECT(eval(GET_SLOT(x, lme4_nlmodelSym), rho));
    if (!isReal(vv) || LENGTH(vv) != n)
	    error(_("evaluated model is not a numeric vector of length %d"), n);
    gg = getAttrib(vv, lme4_gradientSym);
    if (!isReal(gg) || !isMatrix(gg))
	    error(_("gradient attribute of evaluated model must be a numeric matrix"));
    gdims = INTEGER(getAttrib(gg, R_DimSymbol));
    if (gdims[0] != n ||gdims[1] != s)
	    error(_("gradient matrix must be of size %d by %d"), n, s);
    /* colnames of the gradient
     * corresponding to the order of the
     * pnames has been checked */
    Free(eta);
    eta = etaold;
    Memcpy(eta, REAL(vv), n);
    Memcpy(gradient, REAL(gg), ns);
    UNPROTECT(1);
  }
  
  if (muEta) {
    lme4_muEta(mu, muEta, eta, n, dims[lTyp_POS]);
    lme4_varFunc(var, mu, n, dims[vTyp_POS]);
  } else {
    Memcpy(mu, eta, n);
  }
  
  int isLinearModel = !MUETA_SLOT(x) && !V_SLOT(x);
  
  // all of the following are already correct for lmms
  if (isLinearModel) return d[pwrss_POS];
  
  d[wrss_POS] = 0;		/* update resid slot and d[wrss_POS] */
  for (int i = 0; i < n; i++) {
    res[i] = (y[i] - mu[i]) * (srwt ? srwt[i] : 1.0);
    d[wrss_POS] += res[i] * res[i];
  }
  /* store u'u */
  d[usqr_POS] = getSumOfSquares((double*)(cu->x), dims[q_POS]);
  d[pwrss_POS] = d[usqr_POS] + d[wrss_POS];
  // blme edit
  if (canProfileCommonScale(x)) {
    d[sigmaML_POS] = sqrt(d[pwrss_POS]/
                        (srwt ? getSumOfSquares(srwt, n) : (double) n));
    d[sigmaREML_POS] = (gradient || muEta) ? NA_REAL :
      d[sigmaML_POS] * sqrt((((double) n)/((double)(n - p))));
  }
  // blme end
  return d[pwrss_POS];
}

/**
 * Update the contents of the ranef slot in an mer object using the
 * current contents of the u and ST slots.
 *
 * b = T  %*% S %*% t(P) %*% u
 *
 * @param x an mer object
 */
static void update_ranef(SEXP x)
{
  int *Gp = Gp_SLOT(x), *dims = DIMS_SLOT(x), *perm = PERM_VEC(x);
  int nt = dims[nt_POS], q = dims[q_POS];
  double *b = RANEF_SLOT(x), *u = U_SLOT(x), one[] = {1,0};
  int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
  double **st = Alloca(nt, double*);
  R_CheckStack(); 
  
  extractSTParametersAndStructure(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
  /* inverse permutation */
  for (int i = 0; i < q; i++) b[perm[i]] = u[i];
  for (int i = 0; i < nt; i++) {
    for (int k = 0; k < nc[i]; k++) { /* multiply by \tilde{S}_i */
	    double dd = st[i][k * (nc[i] + 1)];
	    int base = Gp[i] + k * nlev[i];
	    for (int kk = 0; kk < nlev[i]; kk++) b[base + kk] *= dd;
    }
    if (nc[i] > 1) {	/* multiply by \tilde{T}_i */
	    F77_CALL(dtrmm)("R", "L", "T", "U", nlev + i, nc + i, one,
                      st[i], nc + i, b + Gp[i], nlev + i);
    }
  }
}

/**
 * We can write the joint distribution of the observations and
 * random effects as:
 *
 *  [ ZLambda  X ]
 *  [ I        0 ]
 *
 * Here we update the cholesky factorizations of the lower-left
 * and lower-right blocks of the above matrix's symmetric cross-
 * product.
 *
 * Update the RZX and RX slots in an mer object. update_L should be
 * called before updateRemainingAugmentedDesignFactors
 *
 * @param regression pointer to an mer object
 */
static void updateRemainingAugmentedDesignMatrixFactors(SEXP regression)
{
  int    *dims      = DIMS_SLOT(regression);
  double *deviances = DEV_SLOT(regression);
  
  double *sqrtModelRowWeight = SXWT_SLOT(regression);
  
  double *denseDesignMatrix         = X_SLOT(regression);
	double *weightedDenseDesignMatrix = (double *) NULL;
  
  int numObservations               = dims[n_POS];
  int numModeledCoefs               = dims[q_POS];
  int numUnmodeledCoefs             = dims[p_POS];
  int numReplicationsOfObservations = dims[s_POS];
  
  double *lowerRightBlockRightFactorization  = RX_SLOT(regression);
  double *offDiagonalBlockRightFactorization = RZX_SLOT(regression);
  
  CHM_SP rotatedSparseDesignMatrix       = A_SLOT(regression);
  CHM_FR upperLeftBlockLeftFactorization = L_SLOT(regression);
  CHM_DN chol_offDiagonalBlockRightFactorization = N_AS_CHM_DN(offDiagonalBlockRightFactorization,
                                                                numModeledCoefs, numUnmodeledCoefs);
  CHM_DN ans;
  R_CheckStack();
  
  // weight X if needed
  if (sqrtModelRowWeight) {			/* Create W^{1/2}GHX in WX */
    weightedDenseDesignMatrix = Calloc(numObservations * numUnmodeledCoefs, double);
    AZERO(weightedDenseDesignMatrix, numObservations * numUnmodeledCoefs);
    
    for (int j = 0; j < numUnmodeledCoefs; ++j) {
	    for (int k = 0; k < numReplicationsOfObservations; ++k) { // will be just 1 for lmm/glmmm
        for (int i = 0; i < numObservations; ++i) {
          int obsIndex    = i + j * numObservations;
          int weightIndex = i + k * numObservations;
          int replIndex   = i + numObservations * (k + j * numReplicationsOfObservations);
          
          weightedDenseDesignMatrix[obsIndex] += sqrtModelRowWeight[weightIndex] * denseDesignMatrix[replIndex];
        }
      }
    }
    
    denseDesignMatrix = weightedDenseDesignMatrix;
    
    /* Replace A by its weighted version, either just the values or entire matrix */
    double *cMatrixExplicitStorage = Cx_SLOT(regression);
    
    if (cMatrixExplicitStorage) rotatedSparseDesignMatrix->x = (void *) cMatrixExplicitStorage;
    else {
	    rotatedSparseDesignMatrix = Cm_SLOT(regression);
	    R_CheckStack();
    }
  }
  
  double d_minusOne[] = { -1.0, 0 }, d_one[] = { 1.0, 0 }, d_zero[] = { 0.0, 0 };
  
  // RZX (temp) = PAX, or
  // temp = Perm * (weighted) rot.Z' * (weighted) X
  multiplyWithPermutation(offDiagonalBlockRightFactorization,
           (int *) upperLeftBlockLeftFactorization->Perm,
           rotatedSparseDesignMatrix,
           denseDesignMatrix,
           numUnmodeledCoefs);
  
  // solve L %*% RZX = PAW^{1/2}GHX, or
  // prodFactor = L^-1 * Perm * rot.Z' * X
  // L being left factor of upper left quadrant of augmented design, rot.Z * rot.Z' + I
  ans = M_cholmod_solve(CHOLMOD_L, upperLeftBlockLeftFactorization,
                        chol_offDiagonalBlockRightFactorization, &cholmodCommon);
  // copy result into the dense, non-cholmod matrix
  Memcpy(offDiagonalBlockRightFactorization, (double *)(ans->x),
         numModeledCoefs * numUnmodeledCoefs);
  M_cholmod_free_dense(&ans, &cholmodCommon);
  
  /* downdate X'X and factor  */
  // compute X'X
  F77_CALL(dsyrk)("U", "T",
                  &numUnmodeledCoefs, &numObservations,
                  d_one,
                  denseDesignMatrix, &numObservations,
                  d_zero,
                  lowerRightBlockRightFactorization, &numUnmodeledCoefs);
  
  // compute X'X - X'X * L^-1 rot.Z'X
  F77_CALL(dsyrk)("U", "T",
                  &numUnmodeledCoefs, &numModeledCoefs,
                  d_minusOne,
                  offDiagonalBlockRightFactorization, &numModeledCoefs,
                  d_one,
                  lowerRightBlockRightFactorization, &numUnmodeledCoefs);

  // take a dense cholesky factorization
  int info;
  F77_CALL(dpotrf)("U",
                   &numUnmodeledCoefs,
                   lowerRightBlockRightFactorization,
                   &numUnmodeledCoefs, &info);
  if (info)
    error(_("Downdated X'X is not positive definite, %d."), info);
    
  // accumulate log(det(RX)^2)
  // matrix is diagonal
  deviances[ldRX2_POS] = 0;
  for (int j = 0; j < numUnmodeledCoefs; ++j) {
    deviances[ldRX2_POS] += 2 * log(lowerRightBlockRightFactorization[j * (numUnmodeledCoefs + 1)]);
  }
  
  if (weightedDenseDesignMatrix) Free(weightedDenseDesignMatrix);
}

/* Level-2 utilities that call at least one of the level-1 utilities */

/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an mer object
 *
 * @return number of iterations to convergence (0 for non-convergence) 
 */
static int update_u(SEXP x)
{
  int *dims = DIMS_SLOT(x);
  int i, n = dims[n_POS], q = dims[q_POS], verb = dims[verb_POS];
  double *Cx = Cx_SLOT(x), *sXwt = SXWT_SLOT(x),
	*res = RESID_SLOT(x), *u = U_SLOT(x),
	cfac = ((double)n) / ((double)q), 
	crit, pwrss, pwrss_old, step;
  double *tmp = Alloca(q, double), *tmp1 = Alloca(q, double),
	*uold = Alloca(q, double), one[] = {1,0}, zero[] = {0,0};
  CHM_FR L = L_SLOT(x);
  CHM_DN cres = N_AS_CHM_DN(res, n, 1),
	ctmp = N_AS_CHM_DN(tmp, q, 1), sol;
  CHM_SP C = Cm_SLOT(x);
  R_CheckStack();
  
  if (!sXwt) return(0);	/* nothing to do for LMMs */
  if (!(L->is_ll)) error(_("L must be LL', not LDL'"));
  //  if (q > n) error(_("q = %d > n = %d"), q, n);
  if (Cx) {		    /* A and C have the same structure */
    C = A_SLOT(x);
    R_CheckStack();
    C->x = (void*)Cx;
  }
  
  /* resetting u to zero at the beginning of each evaluation
   * requires more iterations but is necessary to obtain a
   * repeatable evaluation.  If this is not done the optimization
   * algorithm can take wild steps. */
  AZERO(u, q);
  update_mu(x);
  for (i = 0; ; i++) {
    
    Memcpy(uold, u, q);
    pwrss_old = update_L(x); // update L, W, residuals, A
    /* tmp := PC %*% wtdResid */
    M_cholmod_sdmult(C, 0 /* notrans */, one, zero, cres, ctmp, &cholmodCommon);
    Memcpy(tmp1, tmp, q);
    applyPermutation(tmp, tmp1, (int*)L->Perm, q);
    /* tmp := tmp - u */
    for (int j = 0; j < q; j++) tmp[j] -= u[j];
    /* solve L %*% sol = tmp */
    if (!(sol = M_cholmod_solve(CHOLMOD_L, L, ctmp, &cholmodCommon)))
	    error(_("cholmod_solve (CHOLMOD_L) failed"));
    Memcpy(tmp, (double*)(sol->x), q);
    M_cholmod_free_dense(&sol, &cholmodCommon);
    /* evaluate convergence criterion */
    crit = cfac * getSumOfSquares(tmp, q) / pwrss_old;
    if (crit < CM_TOL) break; /* don't do needless evaluations */
    /* solve t(L) %*% sol = tmp */
    if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &cholmodCommon)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed"));
    Memcpy(tmp, (double*)(sol->x), q);
    M_cholmod_free_dense(&sol, &cholmodCommon);
    
    for (step = 1.0; step > CM_SMIN; step /= 2.0) { /* step halving */
	    for (int j = 0; j < q; j++) u[j] = uold[j] + step * tmp[j];
	    pwrss = update_mu(x);
	    if (verb < 0)
        Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g\n",
                i, step, crit, pwrss, pwrss_old, u[1], u[2]);
	    if (pwrss < pwrss_old) {
        pwrss_old = pwrss;
        break;
	    }
    }
    if (step <= CM_SMIN || i > CM_MAXITER) return 0;
  }
  return i;
}

/**
 * Update the ST, A, L, u, fixef and deviance slots of an mer object.
 *
 * @param x an mer object
 * @param pars double vector of the appropriate length
 *
 */
static void
ST_setPars(SEXP x, const double *pars)
{
  int *Gp = Gp_SLOT(x), nt = DIMS_SLOT(x)[nt_POS], pos = 0;
  int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
  double **st = Alloca(nt, double*);
  R_CheckStack();
  
  extractSTParametersAndStructure(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
  /* install the parameters in the ST slot */
  for (int i = 0; i < nt; i++) {
    int nci = nc[i], ncp1 = nc[i] + 1;
    double *sti = st[i];
    
    for (int j = 0; j < nci; j++)
	    sti[j * ncp1] = pars[pos++];
    for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
        sti[k + j * nci] = pars[pos++];
  }
  
  rotateSparseDesignMatrix(x);
}
/**
 * For a particular model, propagate the changes that result from the
 * selection of new parameters by the optimizer and compute the
 * deviance for that model.
 *
 * @param x an mer object
 * @return updated deviance
 */
static double update_dev(SEXP x)
{
  SEXP flistP = GET_SLOT(x, lme4_flistSym);
  double *d = DEV_SLOT(x), *u = U_SLOT(x);
  int *dims = DIMS_SLOT(x);
  const int q = dims[q_POS], nAGQ = dims[nAGQ_POS], dn = dims[n_POS];
  CHM_FR L = L_SLOT(x);
  
  /* FIXME: This should allow for a GNLMM.  Right now generalized and
   * nonlinear are mutually exclusive.
   */
  /* FIXME: Check these conditions in mer_validate. */
  int isLinearModel = !(MUETA_SLOT(x) || V_SLOT(x));
  
  if (isLinearModel) {
    // shouldn't happen
    error("update_dev called for linear model");
  }
  /* generalized or nonlinear or both */
  update_u(x);
  if (nAGQ < 1) error("nAGQ must be positive");
  if ((nAGQ > 1) & (LENGTH(flistP) != 1))
    error("AGQ method requires a single grouping factor");
  d[ML_POS] = d[ldL2_POS];
  
  if (1/*MUETA_SLOT(x)*/) {	/* GLMM */
    if (nAGQ == 1) {
	    double ans = 0;
	    lme4_devResid(MU_SLOT(x), PWT_SLOT(x), Y_SLOT(x), dims[n_POS],
                    dims[vTyp_POS], &ans, (int*) NULL);
	    d[disc_POS] = ans;
	    d[ML_POS] += d[disc_POS] + d[usqr_POS];
	    return d[ML_POS];
    }
    /* Adaptive Gauss-Hermite quadrature */
    /* Single grouping factor has been checked. */
    const int nl = nlevels(VECTOR_ELT(flistP, 0));
    const int nre = q / nl; 	/* number of random effects per level */
    int *fl0 = INTEGER(VECTOR_ELT(flistP, 0)), *pointer = Alloca(nre, int);
    double *ghw = GHW_SLOT(x), *ghx = GHX_SLOT(x),
    /* store conditional mode */
    *uold = Memcpy(Calloc(q, double), u, q),
    *tmp = Calloc(nl, double),
    w_pro = 1, z_sum = 0;           /* values needed in AGQ evaluation */
    /* constants used in AGQ evaluation */
    const double sigma = (dims[useSc_POS] == 1)?d[sigmaML_POS]:1;
    const double factor = - 0.5 / (sigma * sigma);
    
    R_CheckStack();
    
    AZERO(pointer, nre);
    AZERO(tmp, nl);
    
    while(pointer[nre - 1] < nAGQ){
	    double *z = Calloc(q, double);       /* current abscissas */
	    double *ans = Calloc(nl, double);    /* current penalized residuals in different levels */
      
      /* update abscissas and weights */
	    for(int i = 0; i < nre; ++i){
        for(int j = 0; j < nl; ++j){
          z[i + j * nre] = ghx[pointer[i]];
        }
        w_pro *= ghw[pointer[i]];
        if(!MUETA_SLOT(x))
          z_sum += z[pointer[i]] * z[pointer[i]];
	    }
      
	    CHM_DN cz = N_AS_CHM_DN(z, q, 1), sol;
	    if(!(sol = M_cholmod_solve(CHOLMOD_L, L, cz, &cholmodCommon)))
        error(_("cholmod_solve(CHOLMOD_L) failed"));
	    Memcpy(z, (double *)sol->x, q);
	    M_cholmod_free_dense(&sol, &cholmodCommon);
      
	    for(int i = 0; i < q; ++i) u[i] = uold[i] + sigma * z[i];
	    update_mu(x);
	    
	    AZERO(ans, nl);
	    lme4_devResid(MU_SLOT(x), PWT_SLOT(x), Y_SLOT(x), dims[n_POS],
                    dims[vTyp_POS], ans, fl0);
	    
	    for(int i = 0; i < nre; ++i)
        for(int j = 0; j < nl; ++j)
          ans[j] += u[i + j * nre] * u[i + j * nre];
      
	    for(int i = 0; i < nl; ++i)
        tmp[i] += exp( factor * ans[i] + z_sum) * w_pro / sqrt(PI);
      /* move pointer to next combination of weights and abbsicas */
	    int count = 0;
	    pointer[count]++;
	    while(pointer[count] == nAGQ && count < nre - 1){
        pointer[count] = 0;
        pointer[++count]++;
	    }
      
	    w_pro = 1;
	    z_sum = 0;
      
	    if(z)    Free(z);
	    if(ans)  Free(ans);
	    
    }
    
    for(int j = 0; j < nl; ++j) d[ML_POS] -= 2 * log(tmp[j]);
    Memcpy(u, uold, q);
    update_mu(x);
    d[ML_POS] += MUETA_SLOT(x)?0:(dn * log(2*PI*d[pwrss_POS]/dn));
    if(tmp)   Free(tmp);
    if(uold)  Free(uold);
  } else {  /* NLMM */
    double dn = (double) dims[n_POS];
    
    d[disc_POS] = d[wrss_POS];
    if (nAGQ > 1) {
      /* Adaptive Gauss-Hermite quadrature */
      /* Single grouping factor has been checked. */
	    const int nl = nlevels(VECTOR_ELT(flistP, 0));
	    const int nre = q / nl;
	    int *fl0 = INTEGER(VECTOR_ELT(flistP, 0)), *pointer = Alloca(nre, int);
	    double *ghw = GHW_SLOT(x), *ghx = GHX_SLOT(x),
      *res = RESID_SLOT(x), *tmp = Calloc(nl, double),
      /* store conditional mode */
      *uold = Memcpy(Calloc(q, double), u, q),
      w_pro = 1, z_sum = 0;           /* values needed in AGQ evaluation */
	    const double sigma = d[sigmaML_POS];   /* MLE of sigma */
	    const double factor = - 1 / (2 * sigma * sigma);
      
	    R_CheckStack();
      
	    AZERO(pointer, nre);
	    AZERO(tmp, nl);
	    
	    d[ML_POS] = dn * log(2*PI*d[pwrss_POS]/dn) + d[ldL2_POS];
	    
	    /* implementation of AGQ method (Laplacian will be a trivial case) */
	    AZERO(pointer, nre);                    /* assign initial pointers, all 0 */
	    AZERO(tmp, nl);
      
	    /* add accuracy to integration approximation */
	    while(pointer[nre - 1] < nAGQ){
        double *z = Calloc(q, double);       /* current abscissas */
        double *presid = Calloc(nl, double); /* current penalized residuals in different levels */
        
        /* update abscissas and weights */
        for(int i = 0; i < nre; ++i){
          for(int j = 0; j < nl; ++j){
            z[i + j * nre] = ghx[pointer[i]];
          }
          z_sum += ghx[pointer[i]] * ghx[pointer[i]];
          w_pro *= ghw[pointer[i]];
        }
        CHM_DN cz = N_AS_CHM_DN(z, q, 1), sol;
        if(!(sol = M_cholmod_solve(CHOLMOD_L, L, cz, &cholmodCommon)))
          error(_("cholmod_solve(CHOLMOD_L) failed"));
        Memcpy(z, (double *)sol->x, q);
        M_cholmod_free_dense(&sol, &cholmodCommon);
        for(int i = 0; i < q; ++i){
          u[i] = uold[i] + sigma * z[i];
        }
        update_mu(x);
        
        AZERO(presid, nl);
        for(int i = 0; i < dims[n_POS]; ++i){
          presid[fl0[i]-1] += ( res[i] * res[i] );
        }
        
        for(int i = 0; i < nre; ++i){
          for(int j = 0; j < nl; ++j)
            presid[j] += u[i + j * nre] * u[i + j * nre];
        }
        
        for(int j = 0; j < nl; ++j){
          tmp[j] += exp(factor * presid[j] + z_sum) * w_pro / sqrt(PI);
        }
        
        /* move pointer to next combination of weights and abbsicas */
        int count = 0;
        pointer[count]++;
        while(pointer[count] == nAGQ && count < nre - 1){
          pointer[count] = 0;
          pointer[++count]++;
        }
        if(z) Free(z);
        w_pro = 1;
        z_sum = 0;
        if(presid) Free(presid);
	    }
	    
	    for(int j = 0; j < nl; ++j){
        d[ML_POS] -= ( 2 * log(tmp[j]) );
	    }
	    
	    Memcpy(u, uold, q);
	    update_mu(x);
	    if(tmp)   Free(tmp);
	    if(uold)  Free(uold);
	    
    }
    else{
	    d[ML_POS] = dn*(1 + log(d[pwrss_POS]) + log(2*PI/dn)) + d[ldL2_POS];
    }
  }
  return d[ML_POS];
}

/* Externally callable functions */

/**
 * Create and initialize L
 *
 * @param CmP pointer to the model matrix for the orthogonal random
 * effects (transposed)
 *
 * @return L
 */
SEXP mer_create_L(SEXP CmP)
{
  double one[] = {1, 0};
  CHM_SP Cm = AS_CHM_SP(CmP);
  CHM_FR L;
  R_CheckStack();
  
  L = M_cholmod_analyze(Cm, &cholmodCommon);
  if (!M_cholmod_factorize_p(Cm, one, (int*)NULL, 0 /*fsize*/, L, &cholmodCommon))
    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
          cholmodCommon.status, L->minor, L->n);
  
  return M_chm_factor_to_SEXP(L, 1);
}

/**
 * Generate zeros and weights of Hermite polynomial of order N, for AGQ method
 *
 * changed from fortran in package 'glmmML'
 * @param N order of the Hermite polynomial
 * @param x zeros of the polynomial, abscissas for AGQ
 * @param w weights used in AGQ
 *
 */
static void internal_ghq(int N, double *x, double *w)
{
    int NR, IT, I, K, J;
    double Z = 0, HF = 0, HD = 0;
    double Z0, F0, F1, P, FD, Q, WP, GD, R, R1, R2;
    double HN = 1/(double)N;
    double *X = Calloc(N + 1, double), *W = Calloc(N + 1, double);

    for(NR = 1; NR <= N / 2; NR++){
	if(NR == 1)
	    Z = -1.1611 + 1.46 * sqrt((double)N);
	else
	    Z -= HN * (N/2 + 1 - NR);
	for (IT = 0; IT <= GHQ_MAXIT; IT++) {
	    Z0 = Z;
	    F0 = 1.0;
	    F1 = 2.0 * Z;
	    for(K = 2; K <= N; ++K){
		HF = 2.0 * Z * F1 - 2.0 * (double)(K - 1.0) * F0;
		HD = 2.0 * K * F1;
		F0 = F1;
		F1 = HF;
	    }
	    P = 1.0;
	    for(I = 1; I <= NR-1; ++I){
		P *= (Z - X[I]);
	    }
	    FD = HF / P;
	    Q = 0.0;
	    for(I = 1; I <= NR - 1; ++I){
		WP = 1.0;
		for(J = 1; J <= NR - 1; ++J){
		    if(J != I) WP *= ( Z - X[J] );
		}
		Q += WP;
	    }
	    GD = (HD-Q*FD)/P;
	    Z -= (FD/GD);
	    if (fabs((Z - Z0) / Z) < GHQ_EPS) break;
	}

	X[NR] = Z;
	X[N+1-NR] = -Z;
	R=1.0;
	for(K = 1; K <= N; ++K){
	    R *= (2.0 * (double)K );
	}
	W[N+1-NR] = W[NR] = 3.544907701811 * R / (HD*HD);
    }

    if( N % 2 ){
	R1=1.0;
	R2=1.0;
	for(J = 1; J <= N; ++J){
	    R1=2.0*R1*J;
	    if(J>=(N+1)/2) R2 *= J;
	}
	W[N/2+1]=0.88622692545276*R1/(R2*R2);
	X[N/2+1]=0.0;
    }

    Memcpy(x, X + 1, N);
    Memcpy(w, W + 1, N);

    if(X) Free(X);
    if(W) Free(W);
}

/**
 * Return zeros and weights of Hermite polynomial of order n as a list
 *
 * @param np pointer to a scalar integer SEXP
 * @return a list with two components, the abscissas and the weights.
 *
 */
SEXP lme4_ghq(SEXP np)
{
    int n = asInteger(np);
    SEXP ans = PROTECT(allocVector(VECSXP, 2));

    if (n < 1) n = 1;
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, n ));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, n ));

    internal_ghq(n, REAL(VECTOR_ELT(ans, 0)), REAL(VECTOR_ELT(ans, 1)));
    UNPROTECT(1);
    return ans;
}

/**
 * Book-keeping function exclusively. Not to be exported.
 * 
 * Sets up the state arrays containing options and tolerances
 * for the PORT function S_nlminb_iterate, or non-linear minimization
 * driver.
 *
 * @param regression        an mer object
 * @param stateIV           pointer to array of integers, typically options
 * @param stateIVLength     length of stateIV
 * @param stateV            pointer to array of doubles, typically tolerances
 * @param stateVLength      length of stateV
 * @param optimizationScale pointer to array of doubles, estimate of relative impact of
                            parameter on changes in the target function value
   @param numOptimizationParameters length of optimizationScale
 */
static void initializeOptimizerSettings(SEXP regression,
                                        int *stateIV, int stateIVLength,
                                        double *stateV, int stateVLength,
                                        double *optimizationScale, int numOptimizationParameters)
{
  int *dimensions = DIMS_SLOT(regression);
  
  int isVerbose = dimensions[verb_POS];
  
  if (parametersIncludeUnmodeledCoefs(regression)) {
    // these steps are in lmer.c, but S_Rf_divset undoes these
    // changes so I've left them commented out
    
    // double eta = 1.e-5; // estimated rel. error on computed lpdisc
    
    // stateV[31] = eta;		// RFCTOL
    // stateV[36] = eta;		// SCTOL
    // stateV[41] = eta;		// ETA0
  }
  
  // initialize the state vectors stateV and stateIV
  AZERO(stateIV, stateIVLength);
  AZERO(stateV, stateVLength);
  
  S_Rf_divset(OPT, stateIV, stateIVLength, stateVLength, stateV);

  stateIV[OUTLEV] = (isVerbose < 0) ? -isVerbose : isVerbose; // output 
  stateIV[MXFCAL] = dimensions[mxfn_POS]; // maximum number of function evaluations
  stateIV[MXITER] = dimensions[mxit_POS]; // maximum number of iterations
  
  for (int i = 0; i < numOptimizationParameters; i++) optimizationScale[i] = 1.0;
}

/**
 * Optimize the profiled deviance of an lmer object or the Laplace
 * approximation to the deviance of a nlmer or glmer object.
 *
 * @param regression pointer to an mer object
 *
 * @return R_NilValue
 */
SEXP mer_optimize(SEXP regression)
{
  guaranteeValidPrior(regression); // without a valid prior object, not clear how many parameters to use
  
  int *dims = DIMS_SLOT(regression);
  
  // prepare the stack for optimization
  int numOptimizationParameters = getNumParametersForParameterization(regression, PARAMETERIZATION_PRIOR);
  int numConvergenceParameters  = getNumParametersForParameterization(regression, PARAMETERIZATION_SD_CORRELATION);
  
  // booleans/control
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  int isVerbose = dims[verb_POS];
  int shouldResetCommonScale = FALSE; // !canProfileCommonScale(regression); // should eventually depend on whether or not approximate optimization is in use
  
  // optimization state vectors
  int stateIVLength = S_iv_length(OPT, numOptimizationParameters); // OPT is the optim algorithm type
  int stateVLength  = S_v_length(OPT, numOptimizationParameters);
  
  int    *stateIV = Alloca(stateIVLength, int); // state vectors for optimization routine
  double *stateV  = Alloca(stateVLength, double); // state IV is control, V is for computation
  
  // we optimize on one scale, assess convergence on a second, and let the native lmer code do its
  // thing possibly in a third
  double *optimizationParameters   = Alloca(numOptimizationParameters, double);
  double *convergenceParameters    = Alloca(numConvergenceParameters, double);
  double *oldConvergenceParameters = Alloca(numConvergenceParameters, double);
  
  double *g = (double *) NULL;
  double *h = (double *) NULL;
  double *boxConstraints    = Alloca(2 * numOptimizationParameters, double);
  double *optimizationScale = Alloca(numOptimizationParameters, double); // the parameters are multiplied by this so that they are comparable.
  R_CheckStack();
  
  
  // set initial values
  if (isVerbose) printAllPriors(regression);
  
  initializeOptimizationParameters(regression, optimizationParameters);
  
  initializeOptimizerSettings(regression,
                              stateIV, stateIVLength,
                              stateV,  stateVLength,
                              optimizationScale, numOptimizationParameters);
  
  setBoxConstraints(regression, boxConstraints);
  
  convertOptimizationParametersToConvergence(regression, (const double *) optimizationParameters, convergenceParameters);
  
  double *deviances = DEV_SLOT(regression);
  
  
  // run the optimization
  double commonScaleCache = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  
  double deviance = R_PosInf;
  
  MERCache *cache = (isLinearModel ? createLMMCache(regression) : NULL);
  
  int isMajorIteration;
  int prevIteration = stateIV[NITER];
  int currIteration;
  do {
    currIteration = stateIV[NITER];
    isMajorIteration = (prevIteration != currIteration);
    
    DEBUG_PRINT_ARRAY("par", optimizationParameters, numOptimizationParameters);
        
    if (isMajorIteration) {
      // only check for convergence every major iteration (i.e. not on calls to approximate the gradient)
      prevIteration = currIteration;
      
      Memcpy(oldConvergenceParameters, (const double *) convergenceParameters, numConvergenceParameters);
      convertOptimizationParametersToConvergence(regression, (const double *) optimizationParameters, convergenceParameters);

      if (allApproximatelyEqual(convergenceParameters, oldConvergenceParameters, numConvergenceParameters, 1.0e-10)) {
        stateIV[0] = 6; // absolute function convergence
        break;
      }
    } else if (shouldResetCommonScale) {
      deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS] = commonScaleCache;
    }
    
    if (isAtBoundary(regression, optimizationParameters)) {
      // this call can return true unless a prior with no mass on the boundary is used
      deviance = INFINITY;
    } else {
      // copies the parameters into ST, fixef, and sigma, if necessary
      updateRegressionWithParameters(regression, (const double *) optimizationParameters);

      
      if (isLinearModel) {
        deviance = lmmCalculateDeviance(regression, cache);
      } else {
        rotateSparseDesignMatrix(regression); // this + updateRegression is equal to the old ST_setPars
        deviance = update_dev(regression);
      }
    
      deviance += calculatePriorPenalty(regression, optimizationParameters);
    }
    
    if (isMajorIteration) commonScaleCache = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
        
    S_nlminb_iterate(boxConstraints, optimizationScale, deviance, g, h, stateIV, stateIVLength,
                     stateVLength, numOptimizationParameters, stateV, optimizationParameters);
  } while (stateIV[0] == 1 || stateIV[0] == 2); // continue while no convergence and no error;
  
  DEBUG_PRINT_ARRAY("par", optimizationParameters, numOptimizationParameters);
  
  updateRegressionWithParameters(regression, (const double *) optimizationParameters);
  
  if (isLinearModel) {
    lmmCalculateDeviance(regression, cache);
    deleteLMMCache(cache);
  } else {
    rotateSparseDesignMatrix(regression);
    update_dev(regression);
    // glmms are computed without building other matrix factors as we've
    // optimized over the fixed effects
    updateRemainingAugmentedDesignMatrixFactors(regression);
  }
  
  update_ranef(regression);
  
  dims[cvg_POS] = stateIV[0]; // convergence status
  
  return R_NilValue;
}

/**
 * Extract the parameters from the ST slot of an mer object
 *
 * @param x an mer object
 *
 * @return pointer to a REAL vector
 */
SEXP mer_ST_getPars(SEXP x)
{
  SEXP ans = PROTECT(allocVector(REALSXP, DIMS_SLOT(x)[np_POS]));
  ST_getPars(x, REAL(ans));
  
  UNPROTECT(1); 
  return ans;
}

/**
 * Evaluate starting estimates for the elements of ST
 *
 * @param ST pointers to the nt ST factorizations of the diagonal
 *     elements of Sigma 
 * @param Gpp length nt+1 vector of group pointers for the rows of Zt
 * @param Zt transpose of Z matrix
 *
 */
SEXP mer_ST_initialize(SEXP ST, SEXP Gpp, SEXP Zt)
{
  int *Gp = INTEGER(Gpp),
	*Zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym)),
	*zi = INTEGER(GET_SLOT(Zt, lme4_iSym)), nt = LENGTH(ST);
  int *nc = Alloca(nt, int), *nlev = Alloca(nt, int),
	nnz = INTEGER(GET_SLOT(Zt, lme4_pSym))[Zdims[1]];
  double *rowsqr = Alloca(Zdims[0], double),
	**st = Alloca(nt, double*),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));
  R_CheckStack();
  
  extractSTParametersAndStructure(ST, Gp, st, nc, nlev);
  AZERO(rowsqr, Zdims[0]);
  for (int i = 0; i < nnz; i++) rowsqr[zi[i]] += zx[i] * zx[i];
  for (int i = 0; i < nt; i++) {
    AZERO(st[i], nc[i] * nc[i]);
    for (int j = 0; j < nc[i]; j++) {
	    double *stij = st[i] + j * (nc[i] + 1);
	    for (int k = 0; k < nlev[i]; k++)
        *stij += rowsqr[Gp[i] + j * nlev[i] + k];
	    *stij = sqrt(nlev[i]/(0.375 * *stij));
    }
  }
  return R_NilValue;
}

/**
 * Update the ST slot of an mer object from a REAL vector of
 * parameters and update the sparse model matrix A
 *
 * @param x an mer object
 * @param pars a REAL vector of the appropriate length
 *
 * @return R_NilValue
 */
SEXP mer_ST_setPars(SEXP x, SEXP pars)
{
  SEXP rpars = PROTECT(coerceVector(pars, REALSXP));
  int npar = DIMS_SLOT(x)[np_POS];
  
  if (LENGTH(pars) != npar)
    error(_("pars must be a real vector of length %d"), npar);
  ST_setPars(x, REAL(rpars));
  UNPROTECT(1);
  return R_NilValue;
}

/**
 * Update the whole model for a new set of ST parameters.
 * The only difference between this and a call to update_dev
 * is that the A matrix is also modified.
 *
 * @param regression a bmer object
 *
 * @return updated deviance
 */
SEXP bmer_get_dev(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  
  if (isLinearModel) {
    MERCache *cache = createLMMCache(regression);
    SEXP result = ScalarReal(lmmCalculateDeviance(regression, cache));
    deleteLMMCache(cache);
    return (result);
  }
  
  rotateSparseDesignMatrix(regression);
    
  return(ScalarReal(update_dev(regression)));
}

/**
 * For whatever values of ST and sigma are stored in the model,
 * propagate whatever changes are necessary to compute the deviance.
 * Approximation in this case implies that sigma is not profiled out,
 * and that the approximation is more that less profiling takes place.
 *
 * @param regression a bmer object
 *
 * @return updated deviance
 */
SEXP bmer_approximate_dev(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  
  if (!isLinearModel) error("Approximate deviance not defined for glmm or nlmm.");
  
  MERCache *cache = createLMMCache(regression);
  SEXP result = ScalarReal(lmmApproximateDeviance(regression, cache));
  deleteLMMCache(cache);
  return (result);
}

/**
 * Update the deviance vector in GLMMs, NLMMs and GNLMMs
 * If nAGQ > 1, adaptive Gauss-Hermite quadrature is applied.
 *
 * @param x pointer to an mer object
 *
 * @return R_NilValue
 */
SEXP mer_update_dev(SEXP x)
{
  return ScalarReal(update_dev(x));
}

/**
 * Externally callable update_L.
 * Update the A, L, sqrtrWt and resid slots.  It is assumed that
 * update_mu has already been called at the current values of u and
 * the model parameters.
 *
 * @param x an mer object
 *
 * @return penalized weighted residual sum of squares
 *
 */
SEXP mer_update_L(SEXP x){return ScalarReal(update_L(x));}

/**
 * Externally callable update_mu.
 *
 * Update the eta, v, mu, resid and var slots according to the current
 * values of the parameters and u.  Also evaluate d[wrss_POS] using
 * the current contents of sqrtrWt.  The sqrtrWt slot is updated in
 * update_L.
 *
 * @param x pointer to an mer object
 *
 * @return penalized, weighted residual sum of squares
 */
SEXP mer_update_mu(SEXP x){return ScalarReal(update_mu(x));}

/**
 * Externally callable update_u.
 *
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an mer object
 *
 * @return number of iterations to convergence (0 for non-convergence) 
 */
SEXP mer_update_u(SEXP x){return ScalarInteger(update_u(x));}

/**
 * Externally callable calculateProjection.
 * Create the projections onto the column spaces of the random effects
 * and the fixed effects.
 *
 * @param x an mer object
 *
 * @return a list with two elements, both REAL vectors
 */
SEXP mer_update_projection(SEXP x)
{
  SEXP ans = PROTECT(allocVector(VECSXP, 2));
  int *dims = DIMS_SLOT(x);
  
  SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, dims[q_POS]));
  SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, dims[p_POS]));
  
  calculateProjectionsForSingleArgumentAnova(x, REAL(VECTOR_ELT(ans, 0)), REAL(VECTOR_ELT(ans, 1)));
  
  UNPROTECT(1);
  return ans;
}

/**
 * Externally callable update_ranef.
 * Update the contents of the ranef slot in an mer object.  For a
 * linear mixed model the conditional estimates of the fixed effects
 * and the conditional mode of u are evaluated first.
 *
 * @param x an mer object
 *
 * @return R_NilValue
 */
SEXP mer_update_ranef(SEXP x)
{
  update_ranef(x);
  return R_NilValue;
}

/**
 * Externally callable updateAugmentedDesignFactorization
 *
 * @param x pointer to an mer object
 *
 * @return profiled deviance or REML deviance
 */
SEXP
mer_update_RX(SEXP x)
{
  updateRemainingAugmentedDesignMatrixFactors(x);
  return ScalarReal(DEV_SLOT(x)[ML_POS]);
}

/**
 * Check validity of an mer object
 *
 * @param x Pointer to an mer object
 *
 * @return TRUE if the object is a valid mer object, otherwise a string
 *         that describes the violation.
 */
SEXP mer_validate(SEXP x)
{
  SEXP GpP = GET_SLOT(x, lme4_GpSym),
	ST = GET_SLOT(x, lme4_STSym),
	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym),
	flistP = GET_SLOT(x, lme4_flistSym), asgnP;
  int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP), *asgn;
  const int n = dd[n_POS], nAGQ = dd[nAGQ_POS],
	nt = dd[nt_POS], nfl = LENGTH(flistP),
	p = dd[p_POS], q = dd[q_POS], s = dd[s_POS];
  int nq, nv = n * s;
  CHM_SP Zt = Zt_SLOT(x), A =  A_SLOT(x);
  CHM_FR L = L_SLOT(x);
  char *buf = Alloca(BUF_SIZE + 1, char);
  R_CheckStack();
  /* check lengths */
  if (nt < 1 || LENGTH(ST) != nt)
    return mkString(_("Slot ST must have length dims['nt']"));
  asgnP = getAttrib(flistP, install("assign"));
  if (!isInteger(asgnP) || LENGTH(asgnP) != nt)
    return mkString(_("Slot flist must have integer attribute 'assign' of length dims['nt']"));
  asgn = INTEGER(asgnP);
  if (nAGQ < 1)
    return mkString(_("nAGQ must be positive"));
  if ((nAGQ > 1) & (nfl != 1))
    return mkString(_("AGQ method requires a single grouping factor"));
  
  for (int i = 0; i < nt; i++)
    if (asgn[i] <= 0 || asgn[i] > nfl)
	    return mkString(_("All elements of the assign attribute must be in [1,length(flist)]"));
  if (LENGTH(GpP) != nt + 1)
    return mkString(_("Slot Gp must have length dims['nt'] + 1"));
  
  if (Gp[0] != 0 || Gp[nt] != q)
    return mkString(_("Gp[1] != 0 or Gp[dims['nt'] + 1] != dims['q']"));
  if (LENGTH(devianceP) != (NULLdev_POS + 1) ||
      LENGTH(getAttrib(devianceP, R_NamesSymbol)) != (NULLdev_POS + 1))
    return mkString(_("deviance slot not named or incorrect length"));
  if (LENGTH(dimsP) != (cvg_POS + 1) ||
      LENGTH(getAttrib(dimsP, R_NamesSymbol)) != (cvg_POS + 1))
    return mkString(_("dims slot not named or incorrect length"));
  if (L->n != q || !L->is_ll || !L->is_monotonic)
    return mkString(_("Slot L must be a monotonic LL' factorization of size dims['q']"));
  if (Zt->nrow != q || Zt->ncol != nv)
    return mkString(_("Slot Zt must by dims['q']  by dims['n']*dims['s']"));
  if (A->nrow != q || A->ncol != nv)
    return mkString(_("Slot A must be dims['q']  by dims['n']*dims['s']"));
  if (chkLen(buf, BUF_SIZE, x, lme4_etaSym, n, 0)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_fixefSym, p, 0)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_muEtaSym, n, 1)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_muSym, n, 0)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_offsetSym, n, 1)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_pWtSym, n, 1)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_ranefSym, q, 0)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_residSym, n, 0)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_sqrtrWtSym, n, 1)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_uSym, q, 0)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_VSym, nv, 1)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_varSym, n, 1)) return(mkString(buf));
  if (chkLen(buf, BUF_SIZE, x, lme4_ySym, n, 0)) return(mkString(buf));
  if (chkDims(buf, BUF_SIZE, x, lme4_XSym, nv, p)) return(mkString(buf));
  if (chkDims(buf, BUF_SIZE, x, lme4_RZXSym, q, p)) return(mkString(buf));
  if (chkDims(buf, BUF_SIZE, x, lme4_RXSym, p, p)) return(mkString(buf));
  
  nq = 0;
  for (int i = 0; i < LENGTH(flistP); i++) {
    SEXP fli = VECTOR_ELT(flistP, i);
    if (!isFactor(fli))
	    return mkString(_("flist must be a list of factors"));
    /* 	nq += dm[0] * LENGTH(getAttrib(fli, R_LevelsSymbol)); */
  }
  for (int i = 0; i < nt; i++) {
    SEXP STi = VECTOR_ELT(ST, i);
    int *dm = INTEGER(getAttrib(STi, R_DimSymbol));
    if (!isMatrix(STi) || !isReal(STi) || dm[0] != dm[1])
	    return
      mkString(_("Slot ST must be a list of square numeric matrices"));
    if (Gp[i] > Gp[i + 1])
	    return mkString(_("Gp must be non-decreasing"));
  }
#if 0
  /* FIXME: Need to incorporate the assign attribute in the calculation of nq */
  if (q != nq)
    return mkString(_("q is not sum of columns by levels"));
#endif
  return ScalarLogical(1);
}
