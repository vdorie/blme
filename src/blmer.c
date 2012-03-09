/* This file does:
 *   1) interchange of mer parameterizations between scales useful for
 *     a) evaluating priors
 *     b) assessing convergence
 *     c) the original ST form used by lmer, including:
 *        i) setting initial values for parameters from mer estimate in ST expression
 *       ii) updating ST expression from vector parameters
 *   2) evaluating priors
 *
 * To those ends, most of the code is about going between different parameterizations.
 * All of that should probably be moved elsewhere to streamline the "blmer" bit.
 */
#include <Rmath.h>               /* density functions */

#include "lmer.h"
#include "Syms.h"
#include "blmer.h"
#include "util.h"
#include "wishart.h"
#include "parameters.h"
#include "covariancePrior.h"
#include "unmodeledCoefficientPrior.h"

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

static const char *priorScaleNames[] = {
  "sd", "var", "absolute", "common"
};


// forward declarations
int getNumCovarianceParametersForPrior(priorType_t priorType, int dim);

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

int getNumCovarianceParametersForPrior(priorType_t priorType, int dim)
{
  switch (priorType) {
    case PRIOR_TYPE_CORRELATION:
      return (dim * (dim + 3) / 2);
      break;
    case PRIOR_TYPE_SPECTRAL:
      return (dim * (dim + 1) / 2);
      break;
    case PRIOR_TYPE_DIRECT:
      return (dim * (dim + 1) / 2);
      break;
    case PRIOR_TYPE_NONE:
      return (dim * (dim + 1) / 2);
    default:
      return (0);
      break;
  }
}

void setBoxConstraints(SEXP regression, double *boxConstraints)
{
  // we assume that we optimize on the scale in which we evaluate the prior,
  // as some of our priors are redundantly parameterized
  
  int *dims = DIMS_SLOT(regression);
  SEXP stList = GET_SLOT(regression, lme4_STSym);
  SEXP priorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  int numFactors = dims[nt_POS];
  int numParameters;
  
  for (int i = 0; i < numFactors; ++i) {
    int levelDimension = *INTEGER(getAttrib(VECTOR_ELT(stList, i), R_DimSymbol));
    SEXP prior_i = VECTOR_ELT(priorList, i);
    
    priorType_t priorType = PRIOR_TYPE_SLOT(prior_i);
    
    setCovarianceConstraints(prior_i, boxConstraints, levelDimension);
        
    boxConstraints += 2 * getNumCovarianceParametersForPrior(priorType, levelDimension);
    
  } // level loop end
  
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

void initializeOptimizationParameters(SEXP regression, double *parameters) {
  int *dims = DIMS_SLOT(regression);
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
    
          copySTMatrixToVector(stMatrix, levelDimension, stVector);
          convertSTToPriorCorrelation(stVector, levelDimension, parameters);
        }
        break;
        
      case PRIOR_TYPE_SPECTRAL:
        {
          double *stVector = Alloca(levelDimension * (levelDimension + 1) / 2, double);
          R_CheckStack();
          
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
    double *unmodeledCoefficients = FIXEF_SLOT(regression);
    int numUnmodeledCoefficients = dims[p_POS];
    
    Memcpy(parameters, unmodeledCoefficients, numUnmodeledCoefficients);
    parameters += numUnmodeledCoefficients;
  }
  
  if (parametersIncludeCommonScale(regression)) {
    *parameters++ = DEV_SLOT(regression)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
  } else {
    SEXP commonScalePrior = GET_SLOT(regression, blme_commonScalePriorSym);
    priorType_t priorType = PRIOR_TYPE_SLOT(commonScalePrior);
    
    if (priorType == PRIOR_TYPE_DIRECT &&
        PRIOR_FAMILY_POINT == PRIOR_FAMILIES_SLOT(commonScalePrior)[0]) {
      double commonScale = PRIOR_HYPERPARAMETERS_SLOT(commonScalePrior)[0];
      
      if (PRIOR_SCALES_SLOT(commonScalePrior)[0] == PRIOR_SCALE_VARIANCE) commonScale = sqrt(commonScale);
      
      double* dev = DEV_SLOT(regression);
      dev[sigmaML_POS] = dev[sigmaREML_POS] = commonScale;
    }
  }
}

// this should copy in the parameters from whatever form they're in into
// the ST matrices (and fixef, if necessary). lmer also expects
void updateRegressionWithParameters(SEXP regression, const double *parameters)
{
  int *dims = DIMS_SLOT(regression);
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
// common scale factor
int parametersIncludeUnmodeledCoefs(SEXP regression)
{
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression));
  if (!isLinearModel) return(TRUE);
  
  return(FALSE); // for now, no linear model uses them
}

int parametersIncludeCommonScale(SEXP regression)
{
  return(FALSE);
}

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
    priorScale_t scale = PRIOR_SCALES_SLOT(unmodeledCoefficientPrior)[0];
    
    if (scale == PRIOR_SCALE_COMMON) return(TRUE);
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
      PRIOR_FAMILIES_SLOT(commonScalePrior)[0] != PRIOR_FAMILY_INVGAMMA) {
    // can only profile if is conjugate
    return(FALSE);
  }
  
  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  priorType = PRIOR_TYPE_SLOT(unmodeledCoefficientPrior);
  
  if (priorType == PRIOR_TYPE_NONE) return(FALSE);
  
  if (priorType == PRIOR_TYPE_DIRECT) {
    priorScale_t scale = PRIOR_SCALES_SLOT(unmodeledCoefficientPrior)[0];
    
    if (scale == PRIOR_SCALE_COMMON) return(FALSE);
  }
  
  return(TRUE);
}

/**
 * If there are any priors, they're going to kick back and parameter
 * set that the optimizer chooses where a parameter is at the boundary.
 * This checks against that happening and we simply throw it back to
 * the optimizer.
 */
int isAtBoundary(SEXP regression, double *parameters)
{
  int *dims = DIMS_SLOT(regression);
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

/**
 * For each "level"/term/whatnot, calculate p(Sigma) for the specified
 * type and p(beta) or p(sigma.sq) if necessary.
 */
double calculatePriorPenalty(SEXP regression, double *parameters)
{
  int *dims = DIMS_SLOT(regression);
  SEXP stList = GET_SLOT(regression, lme4_STSym);
  SEXP covariancePriorList = GET_SLOT(regression, blme_covariancePriorSym);
  
  int isLinearModel = !(MUETA_SLOT(regression) || V_SLOT(regression)); 

  double *deviances = DEV_SLOT(regression); // need to use this to get common scale factor (sigma.sq)
  double commonScale = 1.0;
  if (isLinearModel) {
    commonScale = deviances[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS];
    commonScale *= commonScale; // stored as an sd, not a variance
  }
  
  int numFactors = dims[nt_POS];
  
  double deviance = 0.0;
  priorType_t priorType;
  
  for (int i = 0; i < numFactors; ++i) {
    SEXP st_i    = VECTOR_ELT(stList, i);
    SEXP prior_i = VECTOR_ELT(covariancePriorList, i);
    
    priorType = PRIOR_TYPE_SLOT(prior_i);
    int levelDimension = INTEGER(getAttrib(st_i, R_DimSymbol))[0];
    
    deviance   += calculateCovarianceDeviance(prior_i, commonScale, parameters, levelDimension);
    
    parameters += getNumCovarianceParametersForPrior(priorType, levelDimension);
  }
  
  SEXP unmodeledCoefficientPrior = GET_SLOT(regression, blme_unmodeledCoefficientPriorSym);
  priorType = PRIOR_TYPE_SLOT(unmodeledCoefficientPrior);
  
  if (priorType == PRIOR_TYPE_DIRECT &&
      parametersIncludeUnmodeledCoefs(regression))
  {
    int numUnmodeledCoefs = dims[p_POS];
    deviance += calculateUnmodeledCoefficientDeviance(unmodeledCoefficientPrior, commonScale,
                                                      parameters, numUnmodeledCoefs);
    parameters += numUnmodeledCoefs;
  }
  
  SEXP commonScalePrior = GET_SLOT(regression, blme_commonScalePriorSym);
  priorType = PRIOR_TYPE_SLOT(commonScalePrior);
  
  if (priorType == PRIOR_TYPE_DIRECT &&
      parametersIncludeCommonScale(regression))
  {
//    deviance += calclateCommonScaleDeviance(*parameters++, commonScalePrior);
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

SEXP bmer_getScaleEnumeration()
{
  SEXP result = PROTECT(allocVector(STRSXP, PRIOR_SCALE_END - PRIOR_SCALE_SD));
  
  for (priorScale_t scale = PRIOR_SCALE_SD; scale < PRIOR_SCALE_END; ++scale) {
    SET_STRING_ELT(result, (int) scale, mkChar(priorScaleNames[scale]));
  }
  
  UNPROTECT(1);
  
  return (result);
}

SEXP bmer_calculatePriorPenalty(SEXP regression)
{
  int numParameters = getNumParametersForParameterization(regression, PARAMETERIZATION_PRIOR);
  double *parameters = Alloca(numParameters, double);
  R_CheckStack();
  
  initializeOptimizationParameters(regression, parameters);
  
  return(ScalarReal(calculatePriorPenalty(regression, parameters)));
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
  priorScale_t *scales = (priorScale_t *) INTEGER(scalesExpression);
  
  int numScales = LENGTH(scalesExpression);
  if (scales != NULL && numScales > 0) {
    Rprintf("  scales    : %s", priorScaleNames[scales[0]]);
    for (int j = 1; j < numScales; ++j) Rprintf(" %s", priorScaleNames[scales[j]]);
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
