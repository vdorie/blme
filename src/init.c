#include <R_ext/Rdynload.h>
#include <Matrix.h>

#include "lmer.h"
#include "Syms.h" 
#include "sim.h"
#include "blmer.h"
#include "string.h"

#include "test.h"
#include "lmm.h"
#include "unmodeledCoefficientPrior.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(lme4_ghq, 1),
    CALLDEF(mer_ST_getPars, 1),
    CALLDEF(mer_ST_initialize, 3),
    CALLDEF(mer_ST_setPars, 2),
    CALLDEF(mer_update_dev, 1),
    CALLDEF(mer_create_L, 1),
    CALLDEF(mer_optimize, 1),
    CALLDEF(mer_update_ranef, 1),
    CALLDEF(mer_update_mu, 1),
    CALLDEF(mer_validate, 1),
    
    CALLDEF(bmer_sim, 2),
    
    CALLDEF(bmer_getObjectiveFunction, 1),
    CALLDEF(bmer_getPriorPenalty, 1),
  
    CALLDEF(bmer_getObjectiveFunctionForFixedCommonScale, 1),
    CALLDEF(bmer_getCommonScaleDerivatives, 1),
    CALLDEF(bmer_getOptimalCommonScale, 1),
    CALLDEF(bmer_getTypeEnumeration, 0),
    CALLDEF(bmer_getFamilyEnumeration, 0),
    CALLDEF(bmer_getPosteriorScaleEnumeration, 0),
    CALLDEF(bmer_getCommonScaleEnumeration, 0),
    CALLDEF(bmer_getScaleInt, 2),
    
    CALLDEF(bmer_findCommaInList, 1),
    
    CALLDEF(bmer_matrixTest, 0),
    CALLDEF(bmer_parametersTest, 0),
    CALLDEF(bmer_lmmTest, 0),

  {NULL, NULL, 0}
};

/** cholmod_common struct local to the package */
cholmod_common cholmodCommon;

/** Need our own CHOLMOD error handler */
void attribute_hidden
lme4_R_cholmod_error(int status, const char *file, int line, const char *message)
{
    if(status < 0) {
#ifdef Failure_in_matching_Matrix
/* This fails unexpectedly with
 *  function 'cholmod_l_defaults' not provided by package 'Matrix'
 * from ../tests/lmer-1.R 's  (l.68)  lmer(y ~ habitat + (1|habitat*lagoon)
 */
	M_cholmod_defaults(&cholmodCommon);/* <--- restore defaults,
				* as we will not be able to .. */
	cholmodCommon.final_ll = 1;	    /* LL' form of simplicial factorization */
#endif

	error("Cholmod error '%s' at file:%s, line %d", message, file, line);
    }
    else
	warning("Cholmod warning '%s' at file:%s, line %d",
		message, file, line);
}

/** Initializer for blme, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void R_init_blme(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);


    M_R_cholmod_start(&cholmodCommon);
    cholmodCommon.final_ll = 1;	    /* LL' form of simplicial factorization */

    lme4_ASym = install("A");
    lme4_CmSym = install("Cm");
    lme4_CxSym = install("Cx");
    lme4_DimSym = install("Dim");
    lme4_GpSym = install("Gp");
    lme4_LSym = install("L");
    lme4_RXSym = install("RX");
    lme4_RZXSym = install("RZX");
    lme4_STSym = install("ST");
    lme4_VSym = install("V");
    lme4_XSym = install("X");
    lme4_XstSym = install("Xst");
    lme4_ZtSym = install("Zt");
    lme4_devianceSym = install("deviance");
    lme4_dimsSym = install("dims");
    lme4_envSym = install("env");
    lme4_etaSym = install("eta");
    lme4_fixefSym = install("fixef");
    lme4_flistSym = install("flist");
    lme4_ghwSym = install("ghw");
    lme4_ghxSym = install("ghx");
    lme4_gradientSym = install("gradient");
    lme4_iSym = install("i");
    lme4_ncSym = install("nc");
    lme4_nlmodelSym = install("nlmodel");
    lme4_muEtaSym = install("muEta");
    lme4_muSym = install("mu");
    lme4_offsetSym = install("offset");
    lme4_pSym = install("p");
    lme4_permSym = install("perm");
    lme4_pWtSym = install("pWt");
    lme4_ranefSym = install("ranef");
    lme4_residSym = install("resid");
    lme4_sigmaSym = install("sigma");
    lme4_sqrtrWtSym = install("sqrtrWt");
    lme4_sqrtXWtSym = install("sqrtXWt");
    lme4_uSym = install("u");
    lme4_varSym = install("var");
    lme4_xSym = install("x");
    lme4_ySym = install("y");
    blme_covariancePriorSym           = install("cov.prior");
    blme_unmodeledCoefficientPriorSym = install("fixef.prior");
    blme_commonScalePriorSym          = install("var.prior");
    
    blme_prior_typeSym            = install("type");
    blme_prior_familiesSym        = install("families");
    blme_prior_scalesSym          = install("scales");
    blme_prior_hyperparametersSym = install("hyperparameters");
}

/** Finalizer for blmer called upon unloading the package.
 *
 */
void R_unload_blmer(DllInfo *dll){
    M_cholmod_finish(&cholmodCommon);
}
