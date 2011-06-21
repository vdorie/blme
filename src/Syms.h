// Symbols, as they are shared between R and C
// They are defined in Syms.c, and the names are synced
// in init.c

// Also includes a variety of SLOT macros for easy access

#ifndef LME4_SYMS_H
#define LME4_SYMS_H

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Visibility.h>

/** Pointers to symbols initialized in R_init_lme4 */
extern
SEXP attribute_hidden
    lme4_ASym,
    lme4_CmSym,
    lme4_CxSym,
    lme4_DimSym,
    lme4_GpSym,
    lme4_LSym,
    lme4_RXSym,
    lme4_RZXSym,
    lme4_STSym,
    lme4_VSym,
    lme4_XSym,
    lme4_XstSym,
    lme4_ZtSym,
    lme4_devianceSym,
    lme4_dimsSym,
    lme4_envSym,
    lme4_etaSym,
    lme4_fixefSym,
    lme4_flistSym,
    lme4_ghwSym,
    lme4_ghxSym,
    lme4_gradientSym,
    lme4_iSym,
    lme4_muEtaSym,
    lme4_muSym,
    lme4_ncSym,
    lme4_nlmodelSym,
    lme4_offsetSym,
    lme4_pSym,
    lme4_pWtSym,
    lme4_permSym,
    lme4_ranefSym,
    lme4_residSym,
    lme4_sigmaSym,
    lme4_sqrtXWtSym,
    lme4_sqrtrWtSym,
    lme4_uSym,
    lme4_varSym,
    lme4_xSym,
    lme4_ySym,
    blme_covariancePriorSym,
    blme_unmodeledCoefficientPriorSym,
    blme_commonScalePriorSym,
    
    blme_prior_typeSym,
    blme_prior_familiesSym,
    blme_prior_scalesSym,
    blme_prior_hyperparametersSym;


/**
 * Extract the slot named nm from the object obj and return a null pointer
 * if the slot has length 0 or a pointer to the REAL contents.
 *
 * @param obj pointer to an S4 object
 * @param nm pointer to a symbol naming the slot to extract
 * 
 * @return pointer to the REAL contents, if nonzero length, otherwise
 * a NULL pointer 
 *
 */
double *SLOT_REAL_NULL(SEXP obj, SEXP nm);

/** Return the integer pointer to the dims slot */
#define DIMS_SLOT(x) INTEGER(GET_SLOT(x, lme4_dimsSym))

/** Return the double pointer to the deviance slot */
#define DEV_SLOT(x) SLOT_REAL_NULL(x, lme4_devianceSym)

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the A slot and return the pointer. */
#define A_SLOT(x) AS_CHM_SP(GET_SLOT(x, lme4_ASym))

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the Zt slot and return the pointer. */
#define Zt_SLOT(x) AS_CHM_SP(GET_SLOT(x, lme4_ZtSym))

/** Return the integer pointer to the Gp slot */
#define Gp_SLOT(x) INTEGER(GET_SLOT(x, lme4_GpSym))

/** Return the double pointer to the fixef slot */
#define FIXEF_SLOT(x) SLOT_REAL_NULL(x, lme4_fixefSym)

/** Return the double pointer to the muEta slot or (double*) NULL if
 * muEta has length 0) */
#define MUETA_SLOT(x) SLOT_REAL_NULL(x, lme4_muEtaSym)

/** Return the double pointer to the V slot or (double*) NULL if
 * V has length 0) */
#define V_SLOT(x) SLOT_REAL_NULL(x, lme4_VSym)

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the Cm slot and return the pointer. */
#define Cm_SLOT(x) AS_CHM_SP(GET_SLOT(x, lme4_CmSym))

/** Return the double pointer to the Cx slot or (double*) NULL if
 * Cx has length 0) */
#define Cx_SLOT(x) SLOT_REAL_NULL(x, lme4_CxSym)

/** Return the double pointer to the eta slot */
#define ETA_SLOT(x) SLOT_REAL_NULL(x, lme4_etaSym)

/** Return the double pointer to the ghw slot */
#define GHW_SLOT(x) SLOT_REAL_NULL(x, lme4_ghwSym)

/** Return the double pointer to the ghx slot */
#define GHX_SLOT(x) SLOT_REAL_NULL(x, lme4_ghxSym)

/** Allocate (alloca) a cholmod_factor struct, populate it with values
 * from the L slot and return the pointer. */
#define L_SLOT(x) AS_CHM_FR(GET_SLOT(x, lme4_LSym))

/** Return the double pointer to the mu slot */
#define MU_SLOT(x) SLOT_REAL_NULL(x, lme4_muSym)

/** Return the double pointer to the offset slot or (double*) NULL if
 * offset has length 0) */
#define OFFSET_SLOT(x) SLOT_REAL_NULL(x, lme4_offsetSym)

/** Return the integer pointer to the permutation vector in the L slot */
#define PERM_VEC(x) INTEGER(GET_SLOT(GET_SLOT(x, lme4_LSym), lme4_permSym))

/** Return the double pointer to the pWt slot or (double*) NULL if
 * pWt has length 0) */
#define PWT_SLOT(x) SLOT_REAL_NULL(x, lme4_pWtSym)

/** Return the double pointer to the ranef slot or (double*) NULL if
 *  ranef has length 0) */
#define RANEF_SLOT(x) SLOT_REAL_NULL(x, lme4_ranefSym)

/** Residual degrees of freedom */
#define RDF(dims) (dims[n_POS] - (dims[isREML_POS] ? dims[p_POS] : 0))

/** Return the double pointer to the resid slot */
#define RESID_SLOT(x) SLOT_REAL_NULL(x, lme4_residSym)

/** Return the double pointer to the RX slot */
#define RX_SLOT(x) SLOT_REAL_NULL(x, lme4_RXSym)

/** Return the double pointer to the RZX slot */
#define RZX_SLOT(x) SLOT_REAL_NULL(x, lme4_RZXSym)

/** Return the double pointer to the sqrtrWt slot or (double*) NULL if
 *  sqrtrWt has length 0) */
#define SRWT_SLOT(x) SLOT_REAL_NULL(x, lme4_sqrtrWtSym)

/** Return the double pointer to the sqrtXWt slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define SXWT_SLOT(x) SLOT_REAL_NULL(x, lme4_sqrtXWtSym)

/** Return the double pointer to the u slot */
#define U_SLOT(x) SLOT_REAL_NULL(x, lme4_uSym)

/** Return the double pointer to the var slot or (double*) NULL if
 * var has length 0) */
#define VAR_SLOT(x) SLOT_REAL_NULL(x, lme4_varSym)

/** Return the double pointer to the X slot */
#define X_SLOT(x) SLOT_REAL_NULL(x, lme4_XSym)

/** Return the double pointer to the y slot */
#define Y_SLOT(x) SLOT_REAL_NULL(x, lme4_ySym)

#endif // LME4_SYMS_H
