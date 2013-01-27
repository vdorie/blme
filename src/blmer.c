#include "blmer.h"
#include "__merCache.h"

/* This file should be for
 *   general support of lmer's innards
 *   exported functions
 * There's been some bleed through, obviously. Gotta work on that
 */
// #include <Rmath.h>               /* density functions */

#include "lmer.h"
#include "Syms.h"

#include "util.h"
#include "parameters.h"
#include "covariancePrior.h"
#include "unmodeledCoefficientPrior.h"
#include "commonScalePrior.h"
#include "lmer_common.h"
#include "lmm.h"
#include "glmm.h"

#define MIN_PARAMETER    1.0e-10


// these names are exported to R and should match
// enumeration in header
static const char *priorTypeNames[] = {
  "none", "direct", "correlation", "spectral"
};

// these names are exported to R
static const char *priorFamilyNames[] = {
  "flat", "gamma", "inverse.gamma", "wishart", "inverse.wishart",
  "normal", "mvt", "point"
};
static const char *priorFamilyAbbreviations[] = {
  "   1", "gamm", "igam", "wish", "iwsh", "nrml", " mvt", " pnt"
};

static const char *priorPosteriorScaleNames[] = {
  "sd", "var"
};

static const char *priorCommonScaleNames[] = {
  "false", "true"
};

// static const char *priorScaleNames[] = {
//   "sd", "var", "absolute", "common"
// };


// forward declarations
void convertPriorToSDCorrelation(SEXP regression, const double *source, double *target);

void copySTVectorToMatrix(const double *source, int dim, double *target);
void copySTMatrixToVector(const double *source, int dim, double *target);

int getNumParametersForParameterization(SEXP regression, parameterizationType_t type)
{
  int *dims = DIMS_SLOT(regression);
  
  SEXP stList              = GET_SLOT(regression, lme4_STSym);
  SEXP covariancePriorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  int numLevels = dims[nt_POS];
  int numParameters = 0;
  
  for (int i = 0; i < numLevels; ++i) {
    int levelDimension = INTEGER(getAttrib(VECTOR_ELT(stList, i), R_DimSymbol))[0];
    
    SEXP covariancePrior_i = VECTOR_ELT(covariancePriorList, i);
    priorType_t priorType = PRIOR_TYPE_SLOT(covariancePrior_i);
    
    if (type == PARAMETERIZATION_PRIOR) {
      numParameters += getNumCovarianceParametersForPrior(priorType, levelDimension);
    } else {
      numParameters += levelDimension * (levelDimension + 1) / 2;
    }
  }
  
  
  if (parametersIncludeUnmodeledCoefs(regression)) {
    numParameters += dims[p_POS];
  }
  if (parametersIncludeCommonScale(regression)) ++numParameters;
  
  return (numParameters);
}

void convertOptimizationParametersToConvergence(SEXP regression, const double *source, double *target) {
  convertPriorToSDCorrelation(regression, source, target);
}

void setBoxConstraints(SEXP regression, double *boxConstraints)
{
  // we assume that we optimize on the scale in which we evaluate the prior,
  // as some of our priors are redundantly parameterized
  
  const int* dims = DIMS_SLOT(regression);
  SEXP stList = GET_SLOT(regression, lme4_STSym);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  int numFactors = dims[nt_POS];
  int numParameters;
  
  for (int i = 0; i < numFactors; ++i) {
    int levelDimension = INTEGER(getAttrib(VECTOR_ELT(stList, i), R_DimSymbol))[0];
    SEXP prior_i = VECTOR_ELT(priorList, i);
    
    priorType_t priorType = PRIOR_TYPE_SLOT(prior_i);
    
    setCovarianceConstraints(prior_i, boxConstraints, levelDimension);
        
    boxConstraints += 2 * getNumCovarianceParametersForPrior(priorType, levelDimension);
    
  }
  
  if (parametersIncludeUnmodeledCoefs(regression)) {
    numParameters = dims[p_POS];
    for (int j = 0; j < numParameters; ++j) {
      *boxConstraints++ = R_NegInf;
      *boxConstraints++ = R_PosInf;
    }
  }
  if (parametersIncludeCommonScale(regression)) {
    *boxConstraints++ = 0;
    *boxConstraints++ = R_PosInf;
  }
}

void copyParametersFromRegression(SEXP regression, double* parameters) {
  const int* dims = DIMS_SLOT(regression);
  int numFactors = dims[nt_POS];
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  SEXP stList = GET_SLOT(regression, lme4_STSym);
  for (int i = 0; i < numFactors; ++i) {
    SEXP st_i = VECTOR_ELT(stList, i);
    double* stMatrix = REAL(st_i);
    int levelDimension = INTEGER(getAttrib(st_i, R_DimSymbol))[0];
    
    SEXP prior_i = VECTOR_ELT(priorList, i);
    priorType_t priorType = PRIOR_TYPE_SLOT(prior_i);
        
    switch (priorType) {
      case PRIOR_TYPE_CORRELATION:
        {
          double stVector[levelDimension * (levelDimension + 1) / 2];
          
          copySTMatrixToVector(stMatrix, levelDimension, stVector);
          convertSTToPriorCorrelation(stVector, levelDimension, parameters);
        }
        break;
        
      case PRIOR_TYPE_SPECTRAL:
        {
          double stVector[levelDimension * (levelDimension + 1) / 2];
          
          copySTMatrixToVector(stMatrix, levelDimension, stVector);
          convertSTToSpectral(stVector, levelDimension, parameters);
        }
        break;
      
      case PRIOR_TYPE_DIRECT:
        {
          copySTMatrixToVector(stMatrix, levelDimension, parameters);
        }
        break;
      
      default:
        copySTMatrixToVector(stMatrix, levelDimension, parameters);
        break;
    }
    
    parameters += getNumCovarianceParametersForPrior(priorType, levelDimension);
  }
  
  if (parametersIncludeUnmodeledCoefs(regression)) {
    double* unmodeledCoefficients = FIXEF_SLOT(regression);
    int numUnmodeledCoefficients = dims[p_POS];
    
    Memcpy(parameters, unmodeledCoefficients, numUnmodeledCoefficients);
    parameters += numUnmodeledCoefficients;
  }
  
  if (parametersIncludeCommonScale(regression)) {
    *parameters++ = DEV_SLOT(regression)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  } else {
    // what we need to stick somewhere is to insert the common scale for when
    // the prior is fixed to a point. do that here
    SEXP commonScalePrior = GET_SLOT(regression, blme_commonScalePriorSym);
    priorType_t priorType = PRIOR_TYPE_SLOT(commonScalePrior);
    
    if (priorType == PRIOR_TYPE_DIRECT &&
        PRIOR_FAMILY_POINT == PRIOR_FAMILIES_SLOT(commonScalePrior)[0]) {
      double commonScale = PRIOR_HYPERPARAMETERS_SLOT(commonScalePrior)[0];
      
      int priorScale = PRIOR_SCALES_SLOT(commonScalePrior)[0];
      if (getPosteriorScaleBit(priorScale) == PRIOR_POSTERIOR_SCALE_SD) commonScale = sqrt(commonScale);
      
      double* dev = DEV_SLOT(regression);
      dev[sigmaML_POS] = dev[sigmaREML_POS] = commonScale;
    }
  }
}

// this should copy in the parameters from whatever form they're in into
// the ST matrices (and fixef, if necessary). basically the inverse
// of the above
void copyParametersIntoRegression(SEXP regression, const double *parameters)
{
  const int* dims = DIMS_SLOT(regression);
  int numFactors = dims[nt_POS];
  
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  SEXP stList = GET_SLOT(regression, lme4_STSym);
  
  for (int i = 0; i < numFactors; ++i) {
    SEXP st_i = VECTOR_ELT(stList, i);
    double *stMatrix = REAL(st_i);
    int levelDimension = INTEGER(getAttrib(st_i, R_DimSymbol))[0];
    
    SEXP prior_i = VECTOR_ELT(priorList, i);
    priorType_t priorType = PRIOR_TYPE_SLOT(prior_i);
    
    switch (priorType) {
      case PRIOR_TYPE_CORRELATION:
        {
          double *stVector = Alloca(levelDimension * (levelDimension + 1) / 2, double);
          R_CheckStack();
          
          convertPriorCorrelationToST(parameters, levelDimension, stVector);
          copySTVectorToMatrix(stVector, levelDimension, stMatrix);
        }
        break;
        
      case PRIOR_TYPE_SPECTRAL:
        {
          double *stVector = Alloca(levelDimension * (levelDimension + 1) / 2, double);
          R_CheckStack();
          
          convertSpectralToST(parameters, levelDimension, stVector);
          copySTVectorToMatrix(stVector, levelDimension, stMatrix);
        }
        break;
      
      case PRIOR_TYPE_DIRECT:
        {
          copySTVectorToMatrix(parameters, levelDimension, stMatrix);
        }
        break;
      
      default:
        copySTVectorToMatrix(parameters, levelDimension, stMatrix);
        break;
    }
    
    parameters += getNumCovarianceParametersForPrior(priorType, levelDimension);
  }
  
  if (parametersIncludeUnmodeledCoefs(regression)) {
    double *unmodeledCoefficients = FIXEF_SLOT(regression);
    int numUnmodeledCoefficients = dims[p_POS];
    
    Memcpy(unmodeledCoefficients, parameters, numUnmodeledCoefficients);
    parameters += numUnmodeledCoefficients;
  }
  
  if (parametersIncludeCommonScale(regression)) {
    double *deviances = DEV_SLOT(regression);
    if (dims[isREML_POS]) {
      deviances[sigmaREML_POS] = *parameters++;
    } else {
      deviances[sigmaML_POS]   = *parameters++;
    }
  }
}

// st is populated as a vector, but we want it as a matrix
void copySTVectorToMatrix(const double *source, int dim, double *target)
{  
  // copy in diagonal
  for (int i = 0; i < dim; ++i) {
    target[i * (dim + 1)] = *source++;
  }
  
  int offset = 0;
  for (int col = 0; col < dim - 1; ++col) {
    offset = col * (dim + 1) + 1;
    // double scale = target[col * (dim + 1)];
    
    for (int row = col + 1; row < dim; ++row) {
      target[offset++] = *source++;
      //target[col + row * dim] = 0.0;
    }
  }
}

// st lives as a matrix, but we want it is as a vector
void copySTMatrixToVector(const double *source, int dim, double *target)
{
  // copy in diagonal
  for (int i = 0; i < dim; ++i) {
    *target++ = source[i * (dim + 1)];
  }
  
  int offset = 0;
  for (int col = 0; col < dim - 1; ++col) {
    offset = col * (dim + 1) + 1;
    // double scale = source[col * (dim + 1)];
    
    for (int row = col + 1; row < dim; ++row) {
      *target++ = source[offset++]; // / scale;
    }
  }
}

// builds a default with "none" at each level
void guaranteeValidPrior(SEXP regression)
{
  int *dims      = DIMS_SLOT(regression);
  int numLevels = dims[nt_POS];
  
  SEXP covariancePrior = GET_SLOT(regression, blme_covariancePriorSym);
  
  if (covariancePrior == R_NilValue || covariancePrior == NULL) {
    PROTECT(covariancePrior = allocVector(VECSXP, numLevels));
    SET_SLOT(regression, blme_covariancePriorSym, covariancePrior);
    
    for (int i = 0; i < numLevels; ++i) {
      SEXP priorSpecification = PROTECT(priorSpecification = NEW_OBJECT(MAKE_CLASS("bmerPrior")));
      
      SET_SLOT(priorSpecification, blme_prior_typeSym, ScalarInteger(PRIOR_TYPE_NONE));
      
      SET_VECTOR_ELT(covariancePrior, i, priorSpecification);
    }
    UNPROTECT(1 + numLevels);
  } else {
    int unprotectCount = 0;
    int numPriors = LENGTH(covariancePrior);
    
    if (numPriors != numLevels) {
      SEXP oldPrior = covariancePrior;
      
      PROTECT(covariancePrior = allocVector(VECSXP, numLevels));
      SET_SLOT(regression, blme_covariancePriorSym, covariancePrior);
      ++unprotectCount;
      
      for (int i = 0; i < numPriors; ++i) {
        SEXP priorSpecification = VECTOR_ELT(oldPrior, i);
        
        SET_VECTOR_ELT(covariancePrior, i, priorSpecification);
      }
      for (int i = numPriors; i < numLevels; ++i) {
        SET_VECTOR_ELT(covariancePrior, i, R_NilValue);
      }
    }
    
    for (int i = 0; i < numLevels; ++i) {
      SEXP priorSpecification = VECTOR_ELT(covariancePrior, i);
      if (priorSpecification == R_NilValue || priorSpecification == NULL) {
        SEXP priorSpecification = PROTECT(priorSpecification = NEW_OBJECT(MAKE_CLASS("bmerPrior")));
        
        SET_SLOT(priorSpecification, blme_prior_typeSym, ScalarInteger(PRIOR_TYPE_NONE));
        
        SET_VECTOR_ELT(covariancePrior, i, priorSpecification);
        
        ++unprotectCount;
      }
    }
    if (unprotectCount > 0) UNPROTECT(unprotectCount);
  }
  
  
  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  
  if (unmodeledCoefficientPrior == R_NilValue || unmodeledCoefficientPrior == NULL) {
    unmodeledCoefficientPrior = PROTECT(unmodeledCoefficientPrior = NEW_OBJECT(MAKE_CLASS("bmerPrior")));
    SET_SLOT(regression, blme_unmodeledCoefficientPriorSym, unmodeledCoefficientPrior);
    UNPROTECT(1);
  }
  if (LENGTH(GET_SLOT(unmodeledCoefficientPrior, blme_prior_typeSym)) == 0) {
    SET_SLOT(unmodeledCoefficientPrior, blme_prior_typeSym, ScalarInteger(PRIOR_TYPE_NONE));
  } 
  
  
  SEXP commonScalePrior = GET_SLOT(regression, blme_commonScalePriorSym);
  
  if (commonScalePrior == R_NilValue || commonScalePrior == NULL) {
    commonScalePrior = PROTECT(commonScalePrior = NEW_OBJECT(MAKE_CLASS("bmerPrior")));
    SET_SLOT(regression, blme_commonScalePriorSym, commonScalePrior);
    SET_SLOT(commonScalePrior, blme_prior_typeSym, ScalarInteger(PRIOR_TYPE_NONE));
    UNPROTECT(1);
  }
  if (LENGTH(GET_SLOT(commonScalePrior, blme_prior_typeSym)) == 0) {
    SET_SLOT(commonScalePrior, blme_prior_typeSym, ScalarInteger(PRIOR_TYPE_NONE));
  } 
}

// For some models, we have to numerically optimize over the unmodeled
// coefficients. This occurs if a) a glmm, b) the prior over the
// unmodeled coefficients is not a gaussian with variance propto the
// common scale factor (not yet implemented)
int parametersIncludeUnmodeledCoefs(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  if (!isLinearModel) return(TRUE);
  
  return(FALSE); // for now, no linear model uses them
}

// not currently used, will potentially be required for non-conjugate
int parametersIncludeCommonScale(SEXP regression)
{
  return(FALSE);
}

/*
int canProfileCommonScale(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  if (!isLinearModel) return(TRUE); // question doesn't apply if is !lmm
  
  SEXP commonScalePrior = GET_SLOT(regression, blme_commonScalePriorSym);
  priorType_t priorType = PRIOR_TYPE_SLOT(commonScalePrior);
  
  if (priorType == PRIOR_TYPE_DIRECT && 
      PRIOR_FAMILIES_SLOT(commonScalePrior)[0] != PRIOR_FAMILY_INVGAMMA)
  {
    // can only profile if is conjugate
    return FALSE;
  }
  
  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  priorType = PRIOR_TYPE_SLOT(unmodeledCoefficientPrior);
  
  if (priorType == PRIOR_TYPE_NONE) return(TRUE);
  
  if (priorType == PRIOR_TYPE_DIRECT) {
    int scale = PRIOR_SCALES_SLOT(unmodeledCoefficientPrior)[0];
    
    if (PRIOR_SCALE_COMMON(scale) == PRIOR_COMMON_SCALE_TRUE) return(TRUE);
  }
  return(FALSE);
}


int commonScaleRequiresOptimization(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  if (!isLinearModel) return(FALSE); // question doesn't apply if is !lmm
  
  SEXP commonScalePrior = GET_SLOT(regression, blme_commonScalePriorSym);
  priorType_t priorType = PRIOR_TYPE_SLOT(commonScalePrior);
  
  if (priorType == PRIOR_TYPE_DIRECT &&
      PRIOR_FAMILIES_SLOT(commonScalePrior)[0] != PRIOR_FAMILY_INVGAMMA &&
      PRIOR_FAMILIES_SLOT(commonScalePrior)[0] != PRIOR_FAMILY_POINT) {
    // requires optimization if not conjugate and not a point
    return(TRUE);
  }
  
  // at this point, prior, if it exists, is conjugate
  
  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  priorType = PRIOR_TYPE_SLOT(unmodeledCoefficientPrior);
  
  if (priorType == PRIOR_TYPE_NONE) return(FALSE);
  
  if (priorType == PRIOR_TYPE_DIRECT) {
    int scale = PRIOR_SCALES_SLOT(unmodeledCoefficientPrior)[0];
    
    if (PRIOR_SCALE_COMMON(scale) == PRIOR_COMMON_SCALE_TRUE) return(FALSE);
  }
  
  return(TRUE);
} */

/**
 * If there are any priors, they're going to kick back and parameter
 * set that the optimizer chooses where a parameter is at the boundary.
 * This checks against that happening and we simply throw it back to
 * the optimizer.
 *
 * It should definitely be more specific and check that the parameter
 * goes to 0 at 0.
 */
int priorsProhibitParameters(SEXP regression, double *parameters)
{
  const int* dims = DIMS_SLOT(regression);
  SEXP stList    = GET_SLOT(regression, lme4_STSym);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  int numFactors = dims[nt_POS];
  
  for (int i = 0; i < numFactors; ++i) {
    SEXP st_i    = VECTOR_ELT(stList, i);
    SEXP prior_i = VECTOR_ELT(priorList, i);
    
    int levelDimension = INTEGER(getAttrib(st_i, R_DimSymbol))[0];
    priorType_t priorType = PRIOR_TYPE_SLOT(prior_i);
    
    switch (priorType) {
      case PRIOR_TYPE_CORRELATION:
        // checks scales and ST diagonal
        for (int j = 0; j < 2 * levelDimension; ++j) if (parameters[j] <= MIN_PARAMETER) return (TRUE);
        break;
      case PRIOR_TYPE_SPECTRAL:
        // just have to check the eigenvalues
        for (int j = 0; j < levelDimension; ++j) if (parameters[j] <= MIN_PARAMETER) return (TRUE);
        break;
      case PRIOR_TYPE_DIRECT:
        // check that none of the scales are 0
        for (int j = 0; j < levelDimension; ++j) if (parameters[j] <= MIN_PARAMETER) return (TRUE);
        break;
      default:
        break;
    }
    
    parameters += getNumCovarianceParametersForPrior(priorType, levelDimension);
  }
  
  return (FALSE);
}

double getPriorPenalty(SEXP regression, MERCache* cache, const double* parameters)
{
  const int* dims = DIMS_SLOT(regression);
  
  double deviance = cache->priorDevianceConstantPart;
  
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  
  int numParametersUsed = 0;
  deviance   += getCovarianceDevianceVaryingPart(regression, parameters, &numParametersUsed);
  parameters += numParametersUsed;

  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  const double* unmodeledCoefficients = (isLinearModel ? FIXEF_SLOT(regression) : parameters);
  int numUnmodeledCoefficients = dims[p_POS];
  
  if (isLinearModel) {
    double commonScale = DEV_SLOT(regression)[isREML_POS ? sigmaREML_POS : sigmaML_POS];
    commonScale *= commonScale;
    
    deviance += getUnmodeledCoefficientDevianceVaryingPart(unmodeledCoefficientPrior, commonScale, unmodeledCoefficients, numUnmodeledCoefficients);
    
    SEXP commonScalePrior = GET_SLOT(regression, blme_commonScalePriorSym);
    deviance += getCommonScaleDevianceVaryingPart(commonScalePrior, commonScale);
  } else {
    deviance   += getUnmodeledCoefficientDevianceVaryingPart(unmodeledCoefficientPrior, 1.0, unmodeledCoefficients, numUnmodeledCoefficients);
    parameters += numUnmodeledCoefficients;
  }
    
  return (deviance);
}

SEXP bmer_getTypeEnumeration()
{
  SEXP result = PROTECT(allocVector(STRSXP, PRIOR_TYPE_END - PRIOR_TYPE_NONE));
  
  for (priorType_t type = PRIOR_TYPE_NONE; type < PRIOR_TYPE_END; ++type) {
    SET_STRING_ELT(result, (int) type, mkChar(priorTypeNames[type]));
  }
  
  UNPROTECT(1);
  
  return (result);
}


SEXP bmer_getFamilyEnumeration()
{
  SEXP result = PROTECT(allocVector(STRSXP, PRIOR_FAMILY_END - PRIOR_FAMILY_FLAT));
  
  for (priorFamily_t family = PRIOR_FAMILY_FLAT; family < PRIOR_FAMILY_END; ++family) {
    SET_STRING_ELT(result, (int) family, mkChar(priorFamilyNames[family]));
  }
  
  UNPROTECT(1);
  
  return (result);
}

SEXP bmer_getPosteriorScaleEnumeration()
{
  SEXP result = PROTECT(allocVector(STRSXP, PRIOR_POSTERIOR_SCALE_END - PRIOR_POSTERIOR_SCALE_SD));
  
  for (priorPosteriorScale_t scale = PRIOR_POSTERIOR_SCALE_SD; scale < PRIOR_POSTERIOR_SCALE_END; ++scale) {
    SET_STRING_ELT(result, (int) scale, mkChar(priorPosteriorScaleNames[scale]));
  }
  
  UNPROTECT(1);
  
  return (result);
}

SEXP bmer_getCommonScaleEnumeration()
{
  SEXP result = PROTECT(allocVector(STRSXP, PRIOR_COMMON_SCALE_END - PRIOR_COMMON_SCALE_FALSE));
  
  for (priorCommonScale_t scale = PRIOR_COMMON_SCALE_FALSE; scale < PRIOR_COMMON_SCALE_END; ++scale) {
    SET_STRING_ELT(result, (int) scale, mkChar(priorCommonScaleNames[scale]));
  }
  
  UNPROTECT(1);
  
  return (result);
}

SEXP bmer_getScaleInt(SEXP posteriorScale, SEXP commonScale)
{
  int result = 0;
  result |= INTEGER(posteriorScale)[0] * PRIOR_SCALE_POSTERIOR_MASK;
  result |= INTEGER(commonScale)[0]    * PRIOR_SCALE_COMMON_MASK;
  
  return(ScalarInteger(result));
}

SEXP bmer_getPriorPenalty(SEXP regression)
{
  int numParameters = getNumParametersForParameterization(regression, PARAMETERIZATION_PRIOR);
  double* parameters = Alloca(numParameters, double);
  R_CheckStack();
  
  copyParametersFromRegression(regression, parameters);
  
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  MERCache* cache = (isLinearModel ? createLMMCache(regression) : createGLMMCache(regression));
  
  double result = getPriorPenalty(regression, cache, parameters);
  
  if (isLinearModel) deleteLMMCache(cache);
  else deleteGLMMCache(cache);
  
  return(ScalarReal(result));
}


void convertPriorToSDCorrelation(SEXP regression, const double *source, double *target)
{
  int *dims = DIMS_SLOT(regression);
  
  SEXP stList    = GET_SLOT(regression, lme4_STSym);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  int numFactors = dims[nt_POS];
  for (int i = 0; i < numFactors; ++i) {
    SEXP st_i = VECTOR_ELT(stList, i);
    SEXP prior_i = VECTOR_ELT(priorList, i);
    
    int levelDimension = INTEGER(getAttrib(st_i, R_DimSymbol))[0];
    priorType_t priorType = PRIOR_TYPE_SLOT(prior_i);
    
    switch (priorType) {
      case PRIOR_TYPE_CORRELATION:
        convertPriorCorrelationToSDCorrelation(source, levelDimension, target);
        break;
      case PRIOR_TYPE_SPECTRAL:
        convertSpectralToSDCorrelation(source, levelDimension, target);
        break;
      case PRIOR_TYPE_DIRECT:
        target[0] = source[0];
        break;
      default:
        convertCholeskyToSDCorrelation(source, levelDimension, target);
        break;
    }
    
    target += levelDimension * (levelDimension + 1) / 2;
    source += getNumCovarianceParametersForPrior(priorType, levelDimension);
  }
  
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression)); 
  if (!isLinearModel) {
    Memcpy(target, source, dims[p_POS]);
  }
}

static void printType(SEXP prior) {
  priorType_t priorType = PRIOR_TYPE_SLOT(prior);

  Rprintf("  type      : %s\n", priorTypeNames[priorType]);
}

static void printFamilies(SEXP prior) {
  SEXP familiesExpression = GET_SLOT(prior, blme_prior_familiesSym);
  priorFamily_t *families = (priorFamily_t *) INTEGER(familiesExpression);
  
  int numFamilies = LENGTH(familiesExpression);
  if (families != NULL && numFamilies > 0) {
    Rprintf("  families  : %s", priorFamilyAbbreviations[families[0]]);
    for (int j = 1; j < numFamilies; ++j) Rprintf(" %s", priorFamilyAbbreviations[families[j]]);
    Rprintf("\n");
  } else {
    Rprintf("  families  :\n");
  }
}

static void printScales(SEXP prior)
{
  SEXP scalesExpression = GET_SLOT(prior, blme_prior_scalesSym);
  int* scales = INTEGER(scalesExpression);
  
  int numScales = LENGTH(scalesExpression);
  if (scales != NULL && numScales > 0) {
    priorPosteriorScale_t posteriorScale = getPosteriorScaleBit(scales[0]);
    priorCommonScale_t    commonScale    = getCommonScaleBit(scales[0]);
    Rprintf("  scales    : (%s common: %s)", priorPosteriorScaleNames[posteriorScale], (commonScale == PRIOR_COMMON_SCALE_TRUE ? "T" : "F"));
    for (int i = 1; i < numScales; ++i) {
      posteriorScale = getPosteriorScaleBit(scales[i]);
      commonScale    = getCommonScaleBit(scales[i]);
      Rprintf(" (%s common: %s)", priorPosteriorScaleNames[posteriorScale], (commonScale == PRIOR_COMMON_SCALE_TRUE ? "T" : "F"));
    }
    Rprintf("\n");
  } else {
    Rprintf("  scales    :\n");
  }
}

static void printParameters(SEXP prior)
{
  SEXP parametersExpression = GET_SLOT(prior, blme_prior_hyperparametersSym);
  double *parameters = REAL(parametersExpression);
  int numParameters = LENGTH(parametersExpression);
  
  if (parameters != NULL && numParameters > 0) {
    Rprintf("  parameters: %f", parameters[0]);
    for (int j = 1; j < numParameters; ++j) Rprintf(" %f", parameters[j]);
    Rprintf("\n");
  } else {
    Rprintf("  parameters:\n");
  }
}

static void printPrior(SEXP prior)
{
  printType(prior);
  printFamilies(prior);
  printScales(prior);
  printParameters(prior);
}

void printAllPriors(SEXP regression)
{
  int *dims = DIMS_SLOT(regression);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  int numLevels = dims[nt_POS];
  
  Rprintf("covariance prior:\n");
  for (int i = 0; i < numLevels; ++i) {
    SEXP prior_i = VECTOR_ELT(priorList, i);
        
    printPrior(prior_i);
    
    Rprintf("\n");
  }
  
  Rprintf("fixef prior:\n");
  SEXP unmodeledCoefPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  printPrior(unmodeledCoefPrior);
  Rprintf("\n");
  
  Rprintf("var prior:\n");
  SEXP commonScalePrior = GET_SLOT(regression, blme_commonScalePriorSym);
  printPrior(commonScalePrior);
  Rprintf("\n");
}


// Externally callable. Sets up and fills a cache, then returns the computed
// derivatives.
//
// Works on the **deviance** scale, so -2 * logLik.
SEXP bmer_getCommonScaleDerivatives(SEXP regression) {
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  
  if (!isLinearModel) return(ScalarReal(NA_REAL));
  
  MERCache* cache = createLMMCache(regression);
  
  int numParameters = getNumParametersForParameterization(regression, PARAMETERIZATION_PRIOR);
  double parameters[numParameters];
  
  copyParametersFromRegression(regression, parameters);
  
  double firstDerivative;
  double secondDerivative;
  _getCommonScaleDerivatives(regression, cache, parameters, &firstDerivative, &secondDerivative);
  
  deleteLMMCache(cache);
  SEXP resultExp = PROTECT(allocVector(REALSXP, 2));
  double* result = REAL(resultExp);
  result[0] = -2.0 * firstDerivative;
  result[1] = -2.0 * secondDerivative;
  UNPROTECT(1);
  
  return(resultExp);
}

/**
 * Update the whole model for a new set of ST parameters.
 * The only difference between this and a call to update_dev
 * is that the A matrix is also modified. The main idea is that
 * you can plug in your own ST (aka Sigma) matrices and get
 * the objective function back out.
 *
 * @param regression a bmer object
 *
 * @return updated deviance
 */
// extern void printLMMCache(const MERCache* cache);
SEXP bmer_getObjectiveFunction(SEXP regression)
{
  double result;
  
  int numParameters = getNumParametersForParameterization(regression, PARAMETERIZATION_PRIOR);
  double parameters[numParameters];
  
  copyParametersFromRegression(regression, parameters);
  
  rotateSparseDesignMatrix(regression);
  
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  
  if (isLinearModel) {
    MERCache* cache = createLMMCache(regression);
    
    result = lmmGetObjectiveFunction(regression, cache, parameters);
    lmmFinalizeOptimization(regression, cache);
    
    // printLMMCache(cache);
    
    deleteLMMCache(cache);
  } else {
    MERCache* cache = createGLMMCache(regression);
    
    result  = REAL(mer_update_dev(regression))[0];
    result += getPriorPenalty(regression, cache, parameters);
    
    deleteGLMMCache(cache);
  }
  
  return(ScalarReal(result));
}

/**
 * For whatever values of ST and sigma are stored in the model,
 * propagate whatever changes are necessary to compute the deviance.
 *
 * @param regression a bmer object
 *
 * @return updated deviance
 */
SEXP bmer_getObjectiveFunctionForFixedCommonScale(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  
  if (!isLinearModel) error("Common scale not applicable for glmm or nlmm.");
  
  int numParameters = getNumParametersForParameterization(regression, PARAMETERIZATION_PRIOR);
  double parameters[numParameters];
  
  copyParametersFromRegression(regression, parameters);
  
  MERCache* cache = createLMMCache(regression);
  
  SEXP result = ScalarReal(lmmGetObjectiveFunctionForFixedCommonScale(regression, cache, parameters));
  lmmFinalizeOptimization(regression, cache);
  
  deleteLMMCache(cache);
  return (result);
}


/**
 * For whatever values of ST is stored in the model,
 * propagate whatever changes are necessary to compute the maximizer in the
 * common scale.
 *
 * @param regression a bmer object
 *
 * @return the new common scale
 */
SEXP bmer_getOptimalCommonScale(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  
  if (!isLinearModel) error("Common scale not applicable for glmm or nlmm.");
  
  int numParameters = getNumParametersForParameterization(regression, PARAMETERIZATION_PRIOR);
  double parameters[numParameters];
  
  copyParametersFromRegression(regression, parameters);
  
  MERCache* cache = createLMMCache(regression);
  
  SEXP result = ScalarReal(_getOptimalCommonScale(regression, cache, parameters));
  
  deleteLMMCache(cache);
  return (result);
}
