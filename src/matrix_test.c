#include "test.h"
#include "matrix.h"

#include <Rmath.h> // for fabs

#include "util.h" // all approx and all equal functions

#define TEST_TOLERANCE 1.0e-6

static int multiplyMatrices_test();
static int multiplyMatricesWithTranspose_test();
static int singleMatrixCrossproduct_test();
static int singleTriangularMatrixCrossproduct_test();
static int singleMatrixCrossproductWithUpdate_test();
static int applyMatrixToVector_test();
static int applyMatrixToVectorWithUpdate_test();
static int transposeMatrix_test();
static int solveSystem_test();
static int solveSymmetricSystem_test();
static int invertLowerTriangularMatrix_test();
static int invertUnitLowerTriangularMatrix_test();
static int invertUpperTriangularMatrix_test();
static int invertUnitUpperTriangularMatrix_test();
static int getSpectralDecompositionOfPositiveDefiniteMatrix_test();
static int getCholeskyDecomposition_test();
static int getLogDeterminant_test();
static int getLogDeterminantOfTriangularMatrix_test();
static int getLogDeterminantOfPositiveDefiniteMatrix_test();
static int getCayleyTransform_test();
static int getInverseCayleyTransform_test();

// R uses column major format, so these are visually the transposes of the matrices

// used in product and transpose functions
static const double testMatrix1[] = { 0.7947243, 0.326298433, 0.6722109, 0.916388121, 
                                      0.7410866, 0.003083023, 0.3647373, 0.594728339 };
static const int testMatrix1Rows    = 4;
static const int testMatrix1Columns = 2;
static const double testMatrix2[] = { 0.1018805, 0.3691942,
                                      0.6659962, 0.2927602,
                                      0.6686710, 0.9462364 };
static const int testMatrix2Rows    = 2;
static const int testMatrix2Columns = 3;


// used as the RHS matrix in solve
static const double testMatrix3[] = {  0.714434954337776, -0.678789175115526, 0.197810121346265,
                                      -0.166951147839427, -0.142404233571142, 0.32596844155341 };
static const int testMatrix3Rows    = 3;
static const int testMatrix3Columns = 2;

static const double covarianceMatrix[] = { 9.69011322645875, 1.57486048879486, 2.10741986875928,
                                           1.57486048879486, 6.88755872442923, 2.66126660970971,
                                           2.10741986875928, 2.66126660970971, 4.2643587556814 };
static const int covarianceDim = 3;

static const double choleskyMatrix[] = { 3.11289466999106, 0.505915122659575, 0.676996844472522,
                                                      0.0,  2.57519098575107,  0.90042396114631,
                                                      0.0,               0.0,  1.73068504311132 };
static const int choleskyDim = 3;

// runs a suite of unit tests
SEXP bmer_matrixTest() {
  if (!multiplyMatrices_test()) error("multiplyMatrices_test() failed.");
  
  if (!multiplyMatricesWithTranspose_test()) error("multiplyMatricesWithTranspose_test() failed.");
  
  if (!singleMatrixCrossproduct_test()) error("singleMatrixCrossproduct_test() failed.");
  
  if (!singleTriangularMatrixCrossproduct_test()) error("singleTriangularMatrixCrossproduct_test() failed.");
  
  if (!singleMatrixCrossproductWithUpdate_test()) error("singleMatrixCrossproductWithUpdate_test() failed.");
  
  if (!applyMatrixToVector_test()) error("applyMatrixToVector_test() failed.");
  
  if (!applyMatrixToVectorWithUpdate_test()) error("applyMatrixToVectorWithUpdate_test() failed.");
  
  if (!transposeMatrix_test()) error("transposeMatrix_test() failed.");
  
  if (!solveSystem_test()) error("solveSystem_test() failed.");
  
  if (!solveSymmetricSystem_test()) error("solveSymmetricSystem_test() failed.");

  if (!invertLowerTriangularMatrix_test()) error("invertLowerTriangularMatrix_test() failed.");
  
  if (!invertUnitLowerTriangularMatrix_test()) error("invertUnitLowerTriangularMatrix_test() failed.");

  if (!invertUpperTriangularMatrix_test()) error("invertUpperTriangularMatrix_test() failed.");
  
  if (!invertUnitUpperTriangularMatrix_test()) error("invertUnitUpperTriangularMatrix_test() failed.");
  
  if (!getSpectralDecompositionOfPositiveDefiniteMatrix_test()) error("getSpectralDecompositionOfPositiveDefiniteMatrix_test() failed.");
  
  if (!getCholeskyDecomposition_test()) error("getCholeskyDecomposition_test() failed.");
  
  if (!getLogDeterminant_test()) error("getLogDeterminant_test() failed.");
  
  if (!getLogDeterminantOfTriangularMatrix_test()) error("getLogDeterminantOfTriangularMatrix_test() failed.");
  
  if (!getLogDeterminantOfPositiveDefiniteMatrix_test()) error("getLogDeterminantOfPositiveDefiniteMatrix_test() failed.");
  
  if (!getCayleyTransform_test()) error("getCayleyTransform_test() failed.");
  
  if (!getInverseCayleyTransform_test()) error("getInverseCayleyTransform_test() failed.");
  
  return (ScalarLogical(TRUE));
}


static int multiplyMatrices_test() {
  double correctAnswer[] = { 0.35457178346387, 0.0343816817133231, 0.20314407828111, 0.312932333295974,
                             0.74624402508098, 0.218216102874039, 0.55447046989404, 0.784423793782448,
                             1.23265220887754, 0.22110356807718, 0.79461564441162, 1.17551576373053 };

  double result[testMatrix1Rows * testMatrix2Columns];

  multiplyMatrices((double *) testMatrix1, testMatrix1Rows, testMatrix1Columns,
                   (double *) testMatrix2, testMatrix2Rows, testMatrix2Columns,
                   result);

  return (allApproximatelyEqual(result, correctAnswer,
                                testMatrix1Rows * testMatrix2Columns, TEST_TOLERANCE));
}

static int multiplyMatricesWithTranspose_test() {
  double correctAnswer[] = { 0.35457178346387,   0.74624402508098,  1.23265220887754,
                             0.0343816817133231, 0.218216102874039, 0.22110356807718,
                             0.20314407828111,   0.55447046989404,  0.79461564441162,
                             0.312932333295974,  0.784423793782448, 1.17551576373053 };
  
  double result[testMatrix2Columns * testMatrix1Rows];
  
  
  multiplyMatricesWithTranspose((double *) testMatrix2, testMatrix2Rows, testMatrix2Columns, TRUE,
                                (double *) testMatrix1, testMatrix1Rows, testMatrix1Columns, TRUE,
                                result);
                                
  return (allApproximatelyEqual(result, correctAnswer,
                                testMatrix2Columns * testMatrix1Rows, TEST_TOLERANCE));
}

static int singleMatrixCrossproduct_test() {                           
  double correctAnswer1[] = { 2.02969206277747, 1.38014788877641, 1.38014788877641, 1.03595394895137 };
  double result1[testMatrix1Columns * testMatrix1Columns];

  double correctAnswer2[] = { 1.18079606171005, 0.261602080789814, 0.80452426250505, 1.1690211106632,
                              0.261602080789814, 0.106480172409074, 0.220465856800378, 0.300849569050003,
                              0.80452426250505, 0.220465856800378, 0.5849007920901, 0.832925692167063,
                              1.1690211106632, 0.300849569050003, 0.832925692167063, 1.19346898551961 };
  double result2[testMatrix1Rows * testMatrix1Rows];
  
  singleMatrixCrossproduct(testMatrix1, testMatrix1Rows, testMatrix1Columns, result1,
                        FALSE, TRIANGLE_TYPE_BOTH);
  singleMatrixCrossproduct(testMatrix1, testMatrix1Rows, testMatrix1Columns, result2,
                        TRUE, TRIANGLE_TYPE_BOTH);
  
  return (allApproximatelyEqual(result1, correctAnswer1, testMatrix1Columns * testMatrix1Columns, TEST_TOLERANCE) &&
          allApproximatelyEqual(result2, correctAnswer2, testMatrix1Rows * testMatrix1Rows, TEST_TOLERANCE));
}

static int singleTriangularMatrixCrossproduct_test()
{
  const double *lowerTriangle = choleskyMatrix;
  double correctAnswer[] = { 9.69011322645875, 1.57486048879486, 2.10741986875928,
                             1.57486048879486, 6.88755872442922, 2.66126660970971,
                             2.10741986875928, 2.66126660970971, 4.2643587556814 };
  int arrayLength = choleskyDim * choleskyDim;
  double result[arrayLength];
  
  singleTriangularMatrixCrossproduct(lowerTriangle, choleskyDim, TRUE, TRIANGLE_TYPE_LOWER,
                                     result);

  return (allApproximatelyEqual(result, correctAnswer, arrayLength, TEST_TOLERANCE));
}

static int singleMatrixCrossproductWithUpdate_test() {
  double updateMatrix[]  = { 1.45070223605259, 0.867170735190178, 0.867170735190178, 6.18867551188776 };
  double correctAnswer[] = { 3.48039429883006, 2.24731862396659, 2.24731862396659, 7.22462946083913 };
  
  int arrayLength = testMatrix1Columns * testMatrix1Columns;
  
  singleMatrixCrossproductWithUpdate(testMatrix1, testMatrix1Rows, testMatrix1Columns, 1.0,
                                     updateMatrix, FALSE, TRIANGLE_TYPE_BOTH);

  return (allApproximatelyEqual(updateMatrix, correctAnswer, arrayLength, TEST_TOLERANCE));
}

static int applyMatrixToVector_test() {
  int arrayLength = testMatrix1Rows;
  
  double correctAnswer[] = { 0.35457178346387, 0.0343816817133231, 0.20314407828111, 0.312932333295974 };
  double result[arrayLength];
  
  applyMatrixToVector(testMatrix1, testMatrix1Rows, testMatrix1Columns, FALSE,
                      testMatrix2 /* just use first testMatrix1Columns elements of testMatrix2 */,
                      result);

  return (allApproximatelyEqual(result, correctAnswer, arrayLength, TEST_TOLERANCE));
}
  
static int applyMatrixToVectorWithUpdate_test()
{
  int arrayLength = testMatrix1Rows;
  
  double result[] = { 0.490900602890179, 0.346424098825082, 0.138162745162845, 0.249486767454073 };
  double correctAnswer[] = { 0.84547238635405, 0.380805780538405, 0.341306823443955, 0.562419100750047 };
  
  applyMatrixToVectorWithUpdate(testMatrix1, testMatrix1Rows, testMatrix1Columns, FALSE,
                                testMatrix2, 1.0, result);
 
  return (allApproximatelyEqual(result, correctAnswer, arrayLength, TEST_TOLERANCE)); 
}

static int transposeMatrix_test() {
  double correctAnswer[] = { 0.7947243, 0.7410866, 0.326298433, 0.003083023,
                             0.6722109, 0.3647373, 0.916388121, 0.594728339 };
  double result[testMatrix1Rows * testMatrix1Columns];
  
  transposeMatrix(testMatrix1, testMatrix1Rows, testMatrix1Columns, result);
  
  return (allEqual(result, correctAnswer, testMatrix1Columns * testMatrix1Rows));
}

static int solveSystem_test() {
  // a rotation matrix
  double leftHandSide[] = {  0.987576802050028, -0.0130012146128751, -0.156598302900223,
                            -0.0359106006778125, -0.988872097570814, -0.144368983527825,
                            -0.152978720126685, 0.148198998189896, -0.977054025181777 };
  
  
  double correctAnswer[] = {  0.68340774192092, 0.617022240899081, -0.40316039594569,
                             -0.214071757425697, 0.0997551565880077, -0.314052969735081 };

  int arrayLength = testMatrix3Rows * testMatrix3Columns;
  double rightHandSide[arrayLength];
  for (int i = 0; i < arrayLength; ++i) rightHandSide[i] = testMatrix3[i];

  double result[arrayLength];
  
  solveSystem(leftHandSide, testMatrix3Rows, rightHandSide, testMatrix3Columns, result);
  
  return (allApproximatelyEqual(correctAnswer, result, testMatrix3Rows * testMatrix3Columns, TEST_TOLERANCE));
}

static int solveSymmetricSystem_test() {
  int arrayLength = covarianceDim * covarianceDim;
  double leftHandSide[arrayLength];
  for (int i = 0; i < arrayLength; ++i) leftHandSide[i] = covarianceMatrix[i];
  
  arrayLength = testMatrix3Rows * testMatrix3Columns;
  double rightHandSide[arrayLength];
  for (int i = 0; i < arrayLength; ++i) rightHandSide[i] = testMatrix3[i];
  
  
  double correctAnswer[] = {  0.076025442521351,   -0.157264299749081, 0.106959767974372,
                             -0.0359970227480686, -0.0643775124392357, 0.134405907105251  };
  
  
  solveSymmetricSystem(leftHandSide, testMatrix3Rows, rightHandSide, testMatrix3Columns);
  
  return (allApproximatelyEqual(correctAnswer, rightHandSide, arrayLength, TEST_TOLERANCE));
}

static int invertLowerTriangularMatrix_test() {
  int arrayLength = choleskyDim * choleskyDim;
  const double *lowerTriangle = choleskyMatrix;
  double correctAnswer[] = { 0.321244406256403, -0.0631108155061788, -0.092827380400752,
                                           0.0,   0.388320713117262, -0.202031719226983,
                                           0.0,                 0.0,  0.577805883271667 };
  double result[arrayLength];

  invertLowerTriangularMatrix(lowerTriangle, choleskyDim, result);
  
  return (allApproximatelyEqual(result, correctAnswer, arrayLength, TEST_TOLERANCE));
}


static int invertUnitLowerTriangularMatrix_test() {
  int arrayLength = choleskyDim * choleskyDim;
  // R uses column-major matrices, so this is indeed the lower triangle
  double lowerTriangle[] = { 1.0, 0.505915122659575, 0.676996844472522,
                             0.0,               1.0,  0.90042396114631,
                             0.0,               0.0,               1.0 };
  double correctAnswer[] = { 1.0, -0.505915122659575, -0.221458745723566,
                             0.0,                1.0,  -0.90042396114631,
                             0.0,                0.0,                1.0 };
  double result[arrayLength];

  invertUnitLowerTriangularMatrix(lowerTriangle, choleskyDim, result);
  
  return (allApproximatelyEqual(result, correctAnswer, arrayLength, TEST_TOLERANCE));
}

static int invertUpperTriangularMatrix_test() {
  int arrayLength = choleskyDim * choleskyDim;
  double upperTriangle[arrayLength];
  transposeMatrix(choleskyMatrix, choleskyDim, choleskyDim, upperTriangle);
  
  double correctAnswer[] = {   0.321244406256403,                0.0,                0.0,
                             -0.0631108155061788,  0.388320713117262,                0.0,
                              -0.092827380400752, -0.202031719226983,  0.577805883271667 };
  double result[arrayLength];

  invertUpperTriangularMatrix((const double *) upperTriangle, choleskyDim, result);
  
  return (allApproximatelyEqual(result, correctAnswer, arrayLength, TEST_TOLERANCE));
}


static int invertUnitUpperTriangularMatrix_test() {
  int arrayLength = choleskyDim * choleskyDim;
  // R uses column-major matrices, so this is indeed the upper triangle
  double upperTriangle[] = {               1.0,              0.0, 0.0,
                             0.505915122659575,              1.0, 0.0,
                             0.676996844472522, 0.90042396114631, 1.0 };
  double correctAnswer[] = {                1.0,               0.0, 0.0,
                             -0.505915122659575,               1.0, 0.0,
                             -0.221458745723566, -0.90042396114631, 1.0 };
  double result[arrayLength];

  invertUnitUpperTriangularMatrix(upperTriangle, choleskyDim, result);
  
  return (allApproximatelyEqual(result, correctAnswer, arrayLength, TEST_TOLERANCE));
}

static int getSpectralDecompositionOfPositiveDefiniteMatrix_test() {
  int arrayLength = covarianceDim * covarianceDim;
  
  double testMatrix[arrayLength];
  for (int i = 0; i < arrayLength; ++i) testMatrix[i] = covarianceMatrix[i];
  
  double correctVectors[] = {  0.78684109006549,  0.475496108192066,  0.393426676877344,
                              0.598296037951374, -0.744097019704913, -0.297256583170085,
                              0.151403269304993, 0.469279315930783,  -0.869973547691708 };
  double correctValues[] = { 11.6955438642859, 6.68442154276825, 2.4620652995152 };
  
  double eigenVectors[arrayLength];
  double eigenValues[covarianceDim];
  
  int lapackResult = getSpectralDecompositionOfPositiveDefiniteMatrix(testMatrix, covarianceDim, eigenValues, eigenVectors);
  if (lapackResult != 0) {
    if (lapackResult > 0) error("%d eigenvectors failed to converge, check IFAIL array.", lapackResult);
    else error("error in call to LAPACK routine 'dsyevx': argument %d illegal", -lapackResult);
  }
  
  return(allApproximatelyEqual(correctValues, eigenValues, covarianceDim, TEST_TOLERANCE) &&
         allApproximatelyAbsolutelyEqual(correctVectors, eigenVectors, arrayLength, TEST_TOLERANCE));
}

static int getCholeskyDecomposition_test() {
  int arrayLength = covarianceDim * covarianceDim;
  // positive definite matrix
  double testMatrix[arrayLength];
  for (int i = 0; i < arrayLength; ++i) testMatrix[i] = covarianceMatrix[i];
  
  double correctAnswer[] = { 3.11289466999106, 0.505915122659575, 0.676996844472522,
                                            0,  2.57519098575107,  0.90042396114631,
                                            0,                 0,  1.73068504311132 };
  getDenseCholeskyDecomposition(testMatrix, covarianceDim, TRIANGLE_TYPE_LOWER);
  
  for (int col = 1; col < covarianceDim; ++col) {
    for (int row = 0; row < col; ++row) testMatrix[row + col * 3] = 0.0;
  }
    
  return(allApproximatelyEqual(testMatrix, correctAnswer, arrayLength, TEST_TOLERANCE));
}

static int getLogDeterminant_test() {
  double correctAnswer = 5.25998812950527;
  int sign;
  double logDet = getLogDeterminant(covarianceMatrix, covarianceDim, &sign);

  return(fabs(correctAnswer - logDet) <= TEST_TOLERANCE && sign == 1);
}

static int getLogDeterminantOfTriangularMatrix_test() {
  double correctAnswer = 2.62999406475264;
  int sign;
  double logDet = getLogDeterminantOfTriangularMatrix(choleskyMatrix, choleskyDim, &sign);
  
  return(fabs(correctAnswer - logDet) <= TEST_TOLERANCE && sign == 1);
}

static int getLogDeterminantOfPositiveDefiniteMatrix_test() {
  double correctAnswer = 5.25998812950527;
  double logDet = getLogDeterminantOfPositiveDefiniteMatrix(covarianceMatrix, covarianceDim);
  
  return(fabs(correctAnswer - logDet) <= TEST_TOLERANCE);
}

static int getCayleyTransform_test() {
  double testMatrix[] = {                0,  -1.0581370565888,  0.1671810257688,
                           1.0581370565888,                 0, 13.5131086511616,
                          -0.1671810257688, -13.5131086511616,                0 };
  double correctAnswer[] = {  0.987576802050028, -0.0130012146128688, -0.156598302900226,
                             -0.035910600677807,  -0.988872097570814, -0.144368983527825,
                             -0.152978720126689,   0.148198998189895, -0.977054025181777 };

  double result[3 * 3];
  
  getCayleyTransform(testMatrix, 3, result);
  
  return(allApproximatelyEqual(result, correctAnswer, 3 * 3, TEST_TOLERANCE));
}

static int getInverseCayleyTransform_test() {
  // not quite the same as above, but close
  double testMatrix[] = {   0.987576802050028, -0.0130012146128751, -0.156598302900223,
                          -0.0359106006778125,  -0.988872097570814, -0.144368983527825,
                           -0.152978720126685,   0.148198998189896, -0.977054025181777 };
  double correctAnswer[] = {                0,  -1.0581370565888,  0.1671810257688,
                              1.0581370565888,                 0, 13.5131086511616,
                             -0.1671810257688, -13.5131086511616,                0 };

  double result[3 * 3];
  
  getInverseCayleyTransform(testMatrix, 3, result);
  
  return(allApproximatelyEqual(result, correctAnswer, 3 * 3, TEST_TOLERANCE));
}
