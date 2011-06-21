#include <R.h>
#include <Rdefines.h>
#include <R_ext/Visibility.h>

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
double *SLOT_REAL_NULL(SEXP obj, SEXP nm)
{
  SEXP pt = GET_SLOT(obj, nm);
  return LENGTH(pt) ? REAL(pt) : (double*) NULL;
}
