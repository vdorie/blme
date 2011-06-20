#include "test.h"

#include <Matrix.h>

#include "util.h"
#include "Syms.h"
#include "lmer.h"
#include "blmer.h"

// used to create the test regression
static const double testDenseDesignMatrixColumn2[] = { 0.503607972233726, 1.08576936214569, -0.69095383969683, -1.28459935387219, 0.046726172188352, -0.235706556439501, -0.542888255010254, -0.433310317456782, -0.649471646796233, 0.726750747385451, 1.1519117540872, 0.992160365445798, -0.429513109491881, 1.23830410085338, -0.279346281854269, 1.75790308981071, 0.560746090888056, -0.452783972553158, -0.832043296117832, -1.16657054708471, -1.0655905803883, -1.563782051071, 1.15653699715018, 0.83204712857239, -0.227328691424755, 0.266137361672105, -0.376702718583628, 2.44136462889459, -0.795339117255372, -0.0548774737115786, 0.250141322854153, 0.618243293566247, -0.172623502645857, -2.22390027400994, -1.26361438497058, 0.358728895971352, -0.0110454784656636, -0.940649162618608, -0.115825322156954, -0.814968708869917, 0.242263480859686, -1.4250983947325, 0.36594112304922, 0.248412648872596, 0.0652881816716207, 0.0191563916602738, 0.257338377155533, -0.649010077708898, -0.119168762418038, 0.66413569989411 };
static const double testDenseDesignMatrixColumn3[] = { 1.10096910219409, 0.14377148075807, -0.117753598165951, -0.912068366948338, -1.43758624082998, -0.797089525071965, 1.25408310644997, 0.77214218580453, -0.21951562675344, -0.424810283377287, -0.418980099421959, 0.996986860909106, -0.275778029088027, 1.2560188173061, 0.646674390495345, 1.29931230256343, -0.873262111744435, 0.00837095999603331, -0.880871723252545, 0.59625901661066, 0.119717641289537, -0.282173877322451, 1.45598840106634, 0.229019590694692, 0.996543928544126, 0.781859184600258, -0.776776621764597, -0.615989907707918, 0.0465803028049967, -1.13038577760069, 0.576718781896486, -1.28074943178832, 1.62544730346494, -0.500696596002705, 1.67829720781629, -0.412519887482398, -0.97228683550556, 0.0253828675878054, 0.0274753367451927, -1.68018272239593, 1.05375086302862, -1.11959910457218, 0.335617209968815, 0.494795767113158, 0.138052708711737, -0.118792025778828, 0.197684262345795, -1.06869271125479, -0.80321321736474, -1.11376513631953 };

static const int testSparseDesignMatrixNonZeroRowIndices[] = { 4, 9, 12, 17, 22, 1, 6, 13, 18, 23, 3, 8, 10, 15, 20, 3, 8, 11, 16, 21, 4, 9, 13, 18, 23, 4, 9, 10, 15, 20, 3, 8, 10, 15, 20, 2, 7, 13, 18, 23, 4, 9, 12, 17, 22, 4, 9, 12, 17, 22, 1, 6, 12, 17, 22, 2, 7, 12, 17, 22, 3, 8, 10, 15, 20, 4, 9, 13, 18, 23, 1, 6, 12, 17, 22, 1, 6, 12, 17, 22, 4, 9, 11, 16, 21, 0, 5, 12, 17, 22, 4, 9, 12, 17, 22, 1, 6, 11, 16, 21, 4, 9, 13, 18, 23, 0, 5, 12, 17, 22, 0, 5, 13, 18, 23, 0, 5, 12, 17, 22, 0, 5, 12, 17, 22, 0, 5, 10, 15, 20, 1, 6, 12, 17, 22, 0, 5, 13, 18, 23, 2, 7, 12, 17, 22, 4, 9, 12, 17, 22, 1, 6, 12, 17, 22, 0, 5, 11, 16, 21, 3, 8, 12, 17, 22, 4, 9, 12, 17, 22, 3, 8, 11, 16, 21, 1, 6, 12, 17, 22, 1, 6, 12, 17, 22, 1, 6, 12, 17, 22, 3, 8, 13, 18, 23, 2, 7, 10, 15, 20, 2, 7, 13, 18, 23, 2, 7, 12, 17, 22, 0, 5, 12, 17, 22, 4, 9, 12, 17, 22, 0, 5, 13, 18, 23, 3, 8, 14, 19, 24, 1, 6, 11, 16, 21, 2, 7, 10, 15, 20, 0, 5, 12, 17, 22, 1, 6, 12, 17, 22 };
static const int testSparseDesignMatrixIndicesForColumn[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250 };
static const double testSparseDesignMatrixValues[] = { 1, 0.503607972233726, 1, 0.503607972233726, 1.10096910219409, 1, 1.08576936214569, 1, 1.08576936214569, 0.14377148075807, 1, -0.69095383969683, 1, -0.69095383969683, -0.117753598165951, 1, -1.28459935387219, 1, -1.28459935387219, -0.912068366948338, 1, 0.046726172188352, 1, 0.046726172188352, -1.43758624082998, 1, -0.235706556439501, 1, -0.235706556439501, -0.797089525071965, 1, -0.542888255010254, 1, -0.542888255010254, 1.25408310644997, 1, -0.433310317456782, 1, -0.433310317456782, 0.77214218580453, 1, -0.649471646796233, 1, -0.649471646796233, -0.21951562675344, 1, 0.726750747385451, 1, 0.726750747385451, -0.424810283377287, 1, 1.1519117540872, 1, 1.1519117540872, -0.418980099421959, 1, 0.992160365445798, 1, 0.992160365445798, 0.996986860909106, 1, -0.429513109491881, 1, -0.429513109491881, -0.275778029088027, 1, 1.23830410085338, 1, 1.23830410085338, 1.2560188173061, 1, -0.279346281854269, 1, -0.279346281854269, 0.646674390495345, 1, 1.75790308981071, 1, 1.75790308981071, 1.29931230256343, 1, 0.560746090888056, 1, 0.560746090888056, -0.873262111744435, 1, -0.452783972553158, 1, -0.452783972553158, 0.00837095999603331, 1, -0.832043296117832, 1, -0.832043296117832, -0.880871723252545, 1, -1.16657054708471, 1, -1.16657054708471, 0.59625901661066, 1, -1.0655905803883, 1, -1.0655905803883, 0.119717641289537, 1, -1.563782051071, 1, -1.563782051071, -0.282173877322451, 1, 1.15653699715018, 1, 1.15653699715018, 1.45598840106634, 1, 0.83204712857239, 1, 0.83204712857239, 0.229019590694692, 1, -0.227328691424755, 1, -0.227328691424755, 0.996543928544126, 1, 0.266137361672105, 1, 0.266137361672105, 0.781859184600258, 1, -0.376702718583628, 1, -0.376702718583628, -0.776776621764597, 1, 2.44136462889459, 1, 2.44136462889459, -0.615989907707918, 1, -0.795339117255372, 1, -0.795339117255372, 0.0465803028049967, 1, -0.0548774737115786, 1, -0.0548774737115786, -1.13038577760069, 1, 0.250141322854153, 1, 0.250141322854153, 0.576718781896486, 1, 0.618243293566247, 1, 0.618243293566247, -1.28074943178832, 1, -0.172623502645857, 1, -0.172623502645857, 1.62544730346494, 1, -2.22390027400994, 1, -2.22390027400994, -0.500696596002705, 1, -1.26361438497058, 1, -1.26361438497058, 1.67829720781629, 1, 0.358728895971352, 1, 0.358728895971352, -0.412519887482398, 1, -0.0110454784656636, 1, -0.0110454784656636, -0.97228683550556, 1, -0.940649162618608, 1, -0.940649162618608, 0.0253828675878054, 1, -0.115825322156954, 1, -0.115825322156954, 0.0274753367451927, 1, -0.814968708869917, 1, -0.814968708869917, -1.68018272239593, 1, 0.242263480859686, 1, 0.242263480859686, 1.05375086302862, 1, -1.4250983947325, 1, -1.4250983947325, -1.11959910457218, 1, 0.36594112304922, 1, 0.36594112304922, 0.335617209968815, 1, 0.248412648872596, 1, 0.248412648872596, 0.494795767113158, 1, 0.0652881816716207, 1, 0.0652881816716207, 0.138052708711737, 1, 0.0191563916602738, 1, 0.0191563916602738, -0.118792025778828, 1, 0.257338377155533, 1, 0.257338377155533, 0.197684262345795, 1, -0.649010077708898, 1, -0.649010077708898, -1.06869271125479, 1, -0.119168762418038, 1, -0.119168762418038, -0.80321321736474, 1, 0.66413569989411, 1, 0.66413569989411, -1.11376513631953 };

static const double testResponse[] = { 10.2556670398245, 5.77156184409445, 5.2388963179889, -1.01048113754769, -2.40694356326161, 4.87659537730766, 11.0678751360507, 9.95236105245466, -1.21089148312542, 5.80498258814754, 0.684860586400015, 11.0298099678496, 5.71308913175892, 16.7456124781686, 7.4731450675828, 8.51093715973169, 2.82372959757841, 4.3443694315618, -2.96487863376567, 5.88586538048499, 0.976211065859574, 5.62948932886417, 13.6141669016732, 6.07894540274009, 9.57669016288729, 8.56529931086956, 2.02623726516967, 5.335419611075, 5.49337753963837, 0.357086593483832, 5.90421996841909, -0.0078591135781455, 12.5466815499627, -6.66182740887203, 11.9526738840087, 1.5017901567396, 1.15948360785528, 4.97438059345444, 5.91648438771629, 0.262312751711231, 10.5968851648003, -0.281646718883028, 7.38808913061462, 7.91783952393784, 6.24106797266484, 4.86354361177658, 5.04480479030089, 3.28778785631259, 1.36283750127861, -1.0792461964562 };

static const double testObservationWeights[] = { 0.00381387665222907, 0.00253068543212147, 0.0268842109172576, 0.026049490821392, 0.0352712496034507, 0.0117640131347288, 0.0363732801784418, 0.00787598233773136, 0.030897454488213, 0.00875908253475443, 0.00119929989101133, 0.0336135432129502, 0.0267145955150939, 0.0367345999917109, 0.0263537615646393, 0.032876028309309, 0.0141114448422329, 0.0152996366814835, 0.0221359992277507, 0.00371029472391699, 0.00755627846850136, 0.0229306433039551, 0.0293036198379298, 0.0338164933835516, 0.0144975389380483, 0.0311484072530266, 0.00227387001104692, 0.024309809307476, 0.0139066217356676, 0.0229252442673684, 0.0356314693890261, 0.00777690684240461, 0.0143917847174194, 0.0261804195195157, 0.0299524830990098, 0.0203641783916231, 0.0322893691793644, 0.0205532146881175, 0.0195650713155693, 0.0163761413827575, 0.0141272019645207, 0.00481289927179461, 0.0116262989351509, 0.0107885265538658, 0.0308614744202712, 0.0303438490551224, 0.00560673910867697, 0.0201020606053731, 0.0232883713494038, 0.0197244836440227 };

static const double testSTDecompositions[2][9] = { { 1.40987764845173, 1.11529742967833, 0.0, 0.409910474630019 }, { 2.94838215157627, 0.377642085238086, 0.284018269788805, 0, 1.80684813076089, 0.174069288373743, 0, 0, 1.12108512521675 } };

#define TEST_NUM_OBSERVATIONS 50
#define TEST_NUM_UNMODELED_COEFS 3
#define TEST_NUM_FACTORS 2

static const int testNumGroupsPerFactor[] = { 5, 5 };
static const int testNumModeledCoefPerFactor[] = { 2, 3 };

// These values come from an R decomp, namely Cholesky(tcrossprod(C), Imult=1, LDL=FALSE)
// this is not strictly the same as what lmer computes (by avoiding taking the crossproduct it
// arives at a different permutation matrix), so leaving the permutation unchanged the same is necessary
//
// only need the perm and the column counts, as that is all that you'll get from M_cholmod_analyze.
// the rest are so that the SEXP can be converted to a valid CHM_FR
static const int testFactorizationPermutation[] = { 14, 19, 24, 10, 15, 20, 11, 16, 21, 0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 17, 18, 22, 23, 9 };
static const int testFactorizationColumnCounts[] = { 5, 4, 3, 11, 10, 9, 11, 10, 9, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };

static const int testFactorizationNonZeroRowIndices[] = { 0, 1, 2, 12, 16, 1, 2, 12, 16, 2, 12, 16, 3, 4, 5, 9, 11, 12, 13, 15, 16, 17, 24, 4, 5, 9, 11, 12, 13, 15, 16, 17, 24, 5, 9, 11, 12, 13, 15, 16, 17, 24, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 24, 7, 8, 9, 10, 12, 13, 14, 16, 17, 24, 8, 9, 10, 12, 13, 14, 16, 17, 24, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 16, 17, 18, 19, 20, 21, 22, 23, 24, 17, 18, 19, 20, 21, 22, 23, 24, 18, 19, 20, 21, 22, 23, 24, 19, 20, 21, 22, 23, 24, 20, 21, 22, 23, 24, 21, 22, 23, 24, 22, 23, 24, 23, 24, 24 };
static const int testFactorizationIndicesForColumn[] = { 0, 5, 9, 12, 23, 33, 42, 53, 63, 72, 88, 103, 117, 130, 142, 153, 163, 172, 180, 187, 193, 198, 202, 205, 207, 208 };
static const int testFactorizationNumNonZeroes[] = { 5, 4, 3, 11, 10, 9, 11, 10, 9, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
static const int testFactorizationNextColumns[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, -1, 0 };
static const int testFactorizationPrevColumns[] = { 26, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, -1 };
static const int testFactorizationType[] = { CHOLMOD_AMD, TRUE, FALSE, TRUE };

// returns a protected object
SEXP createTestRegression()
{
  SEXP regression = PROTECT(regression = NEW_OBJECT(MAKE_CLASS("bmer")));
  
  int protectCount = 0;
  
  // create and setup the dims slot
  int *dims = INTEGER(ALLOC_SLOT(regression, lme4_dimsSym, INTSXP, (int) (cvg_POS - nt_POS)));
  
  dims[n_POS]  = TEST_NUM_OBSERVATIONS;
  dims[p_POS]  = TEST_NUM_UNMODELED_COEFS;
  dims[nt_POS] = TEST_NUM_FACTORS;
  dims[isREML_POS] = FALSE;
  
  dims[q_POS] = 0;
  for (int i = 0; i < TEST_NUM_FACTORS; ++i) {
    dims[q_POS] += testNumGroupsPerFactor[i] * testNumModeledCoefPerFactor[i];
  }
  dims[np_POS] = dims[q_POS];
  
  int numObservations  = dims[n_POS];
  int numUnmodeledCoef = dims[p_POS];
  int numModeledCoef   = dims[q_POS];
  int numFactors       = dims[nt_POS];
  
  // create the deviance slot
  ALLOC_SLOT(regression, lme4_devianceSym, REALSXP, (int) (NULLdev_POS - ML_POS));
  
  // create and setup the Gp slot
  int *sparseRowsForFactor = INTEGER(ALLOC_SLOT(regression, lme4_GpSym, INTSXP, numFactors + 1));
  
  sparseRowsForFactor[0] = 0;
  for (int i = 0; i < numFactors; ++i) {
    sparseRowsForFactor[i + 1] = testNumGroupsPerFactor[i] * testNumModeledCoefPerFactor[i] + sparseRowsForFactor[i];
  }
  
  // create and setup the X slot
  SEXP denseDesignMatrixExp = ALLOC_SLOT(regression, lme4_XSym, REALSXP, numObservations * numUnmodeledCoef);
  SET_DIMS(denseDesignMatrixExp, numObservations, numUnmodeledCoef);
  double *denseDesignMatrix = REAL(denseDesignMatrixExp);
  for (int i = 0; i < numObservations; ++i) {
    denseDesignMatrix[i]                       = 1.0;
    denseDesignMatrix[i +     numObservations] = testDenseDesignMatrixColumn2[i];
    denseDesignMatrix[i + 2 * numObservations] = testDenseDesignMatrixColumn3[i];
  }
  
  double *response = REAL(ALLOC_SLOT(regression, lme4_ySym, REALSXP, numObservations));
  Memcpy(response, testResponse, numObservations);
  
  // sXwt slot
  double *sqrtObservationWeights = REAL(ALLOC_SLOT(regression, lme4_sqrtXWtSym, REALSXP, numObservations));
  for (int i = 0; i < numObservations; ++i) sqrtObservationWeights[i] = sqrt(testObservationWeights[i]);
  
  // create and setup the Zt slot
  SEXP sparseDesignMatrixExp = PROTECT(sparseDesignMatrixExp = NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
  ++protectCount;
  SET_SLOT(regression, lme4_ZtSym, sparseDesignMatrixExp);
  
  
  int *sdm_dims = INTEGER(ALLOC_SLOT(sparseDesignMatrixExp, install("Dim"), INTSXP, 2));
  sdm_dims[0] = numModeledCoef;
  sdm_dims[1] = numObservations;
  
  int numSparseNonZeroes = 0;
  for (int i = 0; i < numFactors; ++i) numSparseNonZeroes += testNumModeledCoefPerFactor[i];
  numSparseNonZeroes *= numObservations;
  
  int *sdm_nonZeroRowIndices = INTEGER(ALLOC_SLOT(sparseDesignMatrixExp, install("i"), INTSXP, numSparseNonZeroes));
  Memcpy(sdm_nonZeroRowIndices, testSparseDesignMatrixNonZeroRowIndices, numSparseNonZeroes);

  int *sdm_indicesForColumn = INTEGER(ALLOC_SLOT(sparseDesignMatrixExp, install("p"), INTSXP, numObservations + 1));
  Memcpy(sdm_indicesForColumn, testSparseDesignMatrixIndicesForColumn, numObservations + 1);

  double *sdm_values = REAL(ALLOC_SLOT(sparseDesignMatrixExp, install("x"), REALSXP, numSparseNonZeroes));
  Memcpy(sdm_values, testSparseDesignMatrixValues, numSparseNonZeroes);
  
  
  // create and setup the A slot
  SEXP rotatedSparseDesignMatrixExp = PROTECT(rotatedSparseDesignMatrixExp = NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
  ++protectCount;
  SET_SLOT(regression, lme4_ASym, rotatedSparseDesignMatrixExp);
  
  int *rsdm_dims = INTEGER(ALLOC_SLOT(rotatedSparseDesignMatrixExp, install("Dim"), INTSXP, 2));
  rsdm_dims[0] = numModeledCoef;
  rsdm_dims[1] = numObservations;
  
  int *rsdm_nonZeroRowIndices = INTEGER(ALLOC_SLOT(rotatedSparseDesignMatrixExp, install("i"), INTSXP, numSparseNonZeroes));
  Memcpy(rsdm_nonZeroRowIndices, testSparseDesignMatrixNonZeroRowIndices, numSparseNonZeroes);
  
  int *rsdm_indicesForColumn  = INTEGER(ALLOC_SLOT(rotatedSparseDesignMatrixExp, install("p"), INTSXP, numObservations + 1));
  Memcpy(rsdm_indicesForColumn, testSparseDesignMatrixIndicesForColumn, numObservations + 1);
  ALLOC_SLOT(rotatedSparseDesignMatrixExp, install("x"), REALSXP, numSparseNonZeroes);
  
  
  // ST slot
  SEXP stExp = ALLOC_SLOT(regression, lme4_STSym, VECSXP, numFactors);
  for (int i = 0; i < TEST_NUM_FACTORS; ++i) {
    SEXP stExp_i = PROTECT(allocVector(REALSXP, testNumModeledCoefPerFactor[i] * testNumModeledCoefPerFactor[i]));
    ++protectCount;
    SET_VECTOR_ELT(stExp, i, stExp_i);
    SET_DIMS(stExp_i, testNumModeledCoefPerFactor[i], testNumModeledCoefPerFactor[i]);
  
    double *stValues = REAL(stExp_i);
    Memcpy(stValues, testSTDecompositions[i], testNumModeledCoefPerFactor[i] * testNumModeledCoefPerFactor[i]);
  }
  
  // L slot
  SEXP upperLeftBlockLeftFactorizationExp = PROTECT(NEW_OBJECT(MAKE_CLASS("dCHMsimpl")));
  ++protectCount;
  SET_SLOT(regression, lme4_LSym, upperLeftBlockLeftFactorizationExp);
  
  int *ulfblf_permutation = INTEGER(ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("perm"), INTSXP, numModeledCoef));
  Memcpy(ulfblf_permutation, testFactorizationPermutation, numModeledCoef);

  int *ulfblf_columnCounts = INTEGER(ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("colcount"), INTSXP, numModeledCoef));
  Memcpy(ulfblf_columnCounts, testFactorizationColumnCounts, numModeledCoef);
  
  int numFactorizationNonZeroes = 0;
  for (int i = 0; i < numModeledCoef; ++i) numFactorizationNonZeroes += ulfblf_columnCounts[i];
  
  ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("x"), REALSXP, numFactorizationNonZeroes);
  
  int *ulfblf_indicesForColumn = INTEGER(ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("p"), INTSXP, numModeledCoef + 1));
  Memcpy(ulfblf_indicesForColumn, testFactorizationIndicesForColumn, numModeledCoef + 1);
  
  int *ulfblf_nonZeroRowIndices = INTEGER(ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("i"), INTSXP, numFactorizationNonZeroes));
  Memcpy(ulfblf_nonZeroRowIndices, testFactorizationNonZeroRowIndices, numFactorizationNonZeroes);
  
  int *ulfblf_numNonZeroes = INTEGER(ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("nz"), INTSXP, numModeledCoef));
  Memcpy(ulfblf_numNonZeroes, testFactorizationNumNonZeroes, numModeledCoef);
  
  int *ulfblf_nextColumns = INTEGER(ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("nxt"), INTSXP, numModeledCoef + 2));
  Memcpy(ulfblf_nextColumns, testFactorizationNextColumns, numModeledCoef + 2);
  
  int *ulfblf_prevColumns = INTEGER(ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("prv"), INTSXP, numModeledCoef + 2));
  Memcpy(ulfblf_prevColumns, testFactorizationPrevColumns, numModeledCoef + 2);
  
  int *ulfblf_type = INTEGER(ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("type"), INTSXP, 4));
  Memcpy(ulfblf_type, testFactorizationType, 4);
  
  int *ulfblf_dims = INTEGER(ALLOC_SLOT(upperLeftBlockLeftFactorizationExp, install("Dim"), INTSXP, 2));
  ulfblf_dims[0] = ulfblf_dims[1] = numModeledCoef;
  
  // misc slots
  ALLOC_SLOT(regression, lme4_offsetSym, REALSXP, 0);
  ALLOC_SLOT(regression, lme4_varSym,    REALSXP, 0);
  ALLOC_SLOT(regression, lme4_fixefSym,  REALSXP, numUnmodeledCoef);
  ALLOC_SLOT(regression, lme4_uSym,      REALSXP, numModeledCoef);
  ALLOC_SLOT(regression, lme4_CxSym,     REALSXP, numSparseNonZeroes);
  
  SEXP offDiagonalBlockRightFactorizationExp =
    ALLOC_SLOT(regression, lme4_RXSym, REALSXP, numUnmodeledCoef * numUnmodeledCoef);
  AZERO(REAL(offDiagonalBlockRightFactorizationExp), numUnmodeledCoef * numUnmodeledCoef);
  SET_DIMS(offDiagonalBlockRightFactorizationExp, numUnmodeledCoef, numUnmodeledCoef);

  SEXP lowerRightBlockRightFactorizationExp = 
    ALLOC_SLOT(regression, lme4_RZXSym, REALSXP, numModeledCoef * numUnmodeledCoef);
  SET_DIMS(lowerRightBlockRightFactorizationExp, numModeledCoef, numUnmodeledCoef);
  
  guaranteeValidPrior(regression);
  
  // at this point, everything should be jammed into the regression
  // or its objects
  UNPROTECT(protectCount);
  
  return (regression);
}
