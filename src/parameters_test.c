#include "test.h"
#include "parameters.h"
#include "util.h" // all approx and all equal functions

#define TEST_TOLERANCE 1.0e-6

static int convertSDCorrelationToCovariance_test();
static int convertSpectralToCovariance_test();
static int convertSTToCovariance_test();
static int convertSTToCovarianceWithScale_test();
static int convertCholeskyToCovariance_test();
static int convertPriorCorrelationToCovariance_test();
static int convertSTToCovarianceInverse_test();
static int convertSTToCovarianceInverseWithScale_test();

static int convertCovarianceToSDCorrelation_test();
static int convertSpectralToSDCorrelation_test();
static int convertSTToSDCorrelation_test();
static int convertCholeskyToSDCorrelation_test();
static int convertPriorCorrelationToSDCorrelation_test();

static int convertCovarianceToSpectral_test();
static int convertSDCorrelationToSpectral_test();
static int convertSTToSpectral_test();
static int convertCholeskyToSpectral_test();
static int convertPriorCorrelationToSpectral_test();

static int convertCovarianceToST_test();
static int convertSDCorrelationToST_test();
static int convertSpectralToST_test();
static int convertCholeskyToST_test();
static int convertPriorCorrelationToST_test();

static int convertCovarianceToCholesky_test();
static int convertSDCorrelationToCholesky_test();
static int convertSpectralCholesky_test();
static int convertSTToCholesky_test();
static int convertPriorCorrelationToCholesky_test();

static int convertCovarianceToPriorCorrelation_test();
static int convertSDCorrelationToPriorCorrelation_test();
static int convertSpectralToPriorCorrelation_test();
static int convertSTToPriorCorrelation_test();
static int convertCholeskyToPriorCorrelation_test();

// R uses column major format, so these are visually the transposes of the matrices

static const double covarianceMatrix[] = {    8.3427288108994, -0.58085137646235, -0.131753109528372,
                                            -0.58085137646235,  1.19615454427737,  -1.13268851208116,
                                           -0.131753109528372, -1.13268851208116,   7.87952065974224  };
static const int covarianceDim = 3;
static const double sdCorrelationArray[] = {  2.88837823196676,   1.09368850422658,   2.80704838927694,
                                            -0.183872724959694, -0.016250132679827, -0.368949519687758  };
static const int sdCorrelationDim = 3;

static const double spectralArray[] = {  8.39443601351161,   8.0624530711771, 0.961514930230293,    // eigen values
                                        0.101413800853997, 0.023216412427196, -1.18403135576557  }; // skew-symmetric
static const int spectralDim = 3;

static const double stArray[] = {    2.88837823196676,    1.07504118116342,  2.59793459589513,   // S
                                  -0.0696236674627963, -0.0157925676975431, -0.98801442231373 }; // T
static const int stDim = 3;

static const double choleskyMatrix[] = { 2.88837823196676, -0.201099485529233, -0.0456149087644448,
                                                        0,   1.07504118116342,   -1.06215619157065,
                                                        0,                  0,    2.59793459589513  };
static const int choleskyDim = 3;

static const double priorCorrelationArray[] = {   1.69952294246555,    1.03684192679667,   1.61181096779217,    // D
                                                  1.69952294246555,    1.03684192679667,   1.61181096779217,    // S
                                                -0.114122526427135, -0.0166519720108243, -0.635567568284938  }; // T
static const int priorCorrelationDim = 3;

static const double covarianceMatrixInverse[] = {   0.12511917535566, 0.0726248077043794, 0.0125319965588894,
                                                  0.0726248077043795,      1.00990023249,  0.146388465205345,
                                                  0.0125319965588894,  0.146388465205345,  0.148164299932518  };
static const int covarianceInverseDim = 3;

SEXP bmer_parametersTest()
{
  int result = convertSDCorrelationToCovariance_test();
  if (!result) error("convertSDCorrelationToCovariance_test() failed.");
  
  result = convertSpectralToCovariance_test();
  if (!result) error("convertSpectralToCovariance_test() failed.");
  
  result = convertSTToCovariance_test();
  if (!result) error("convertSTToCovariance_test() failed.");
  
  result = convertSTToCovarianceWithScale_test();
  if (!result) error("convertSTToCovarianceWithScale_test() failed.");
  
  result = convertCholeskyToCovariance_test();
  if (!result) error("convertCholeskyToCovariance_test() failed.");
  
  result = convertPriorCorrelationToCovariance_test();
  if (!result) error("convertPriorCorrelationToCovariance_test() failed.");
  
  result = convertSTToCovarianceInverse_test();
  if (!result) error("convertSTToCovarianceInverse_test() failed.");
  
  result = convertSTToCovarianceInverseWithScale_test();
  if (!result) error("convertSTToCovarianceInverseWithScale_test() failed.");
  
  result = convertCovarianceToSDCorrelation_test();
  if (!result) error("convertCovarianceToSDCorrelation_test() failed.");
  
  result = convertSpectralToSDCorrelation_test();
  if (!result) error("convertSpectralToSDCorrelation_test() failed.");
  
  result = convertSTToSDCorrelation_test();
  if (!result) error("convertSTToSDCorrelation_test() failed.");
  
  result = convertCholeskyToSDCorrelation_test();
  if (!result) error("convertCholeskyToSDCorrelation_test() failed.");
  
  result = convertPriorCorrelationToSDCorrelation_test();
  if (!result) error("convertPriorCorrelationToSDCorrelation_test() failed.");
  
  result = convertCovarianceToSpectral_test();
  if (!result) error("convertCovarianceToSpectral_test() failed.");
  
  result = convertSDCorrelationToSpectral_test();
  if (!result) error("convertSDCorrelationToSpectral_test() failed.");
  
  result = convertSTToSpectral_test();
  if (!result) error("convertSTToSpectral_test() failed.");
  
  result = convertCholeskyToSpectral_test();
  if (!result) error("convertCholeskyToSpectral_test() failed.");
  
  result = convertPriorCorrelationToSpectral_test();
  if (!result) error("convertPriorCorrelationToSpectral_test() failed.");
  
  result = convertCovarianceToST_test();
  if (!result) error("convertCovarianceToST_test() failed.");
  
  result = convertSDCorrelationToST_test();
  if (!result) error("convertSDCorrelationToST_test() failed.");
  
  result = convertSpectralToST_test();
  if (!result) error("convertSpectralToST_test() failed.");
  
  result = convertCholeskyToST_test();
  if (!result) error("convertCholeskyToST_test() failed.");
  
  result = convertPriorCorrelationToST_test();
  if (!result) error("convertPriorCorrelationToST_test() failed.");
  
  result = convertCovarianceToCholesky_test();
  if (!result) error("convertCovarianceToCholesky_test() failed.");
  
  result = convertSDCorrelationToCholesky_test();
  if (!result) error("convertSDCorrelationToCholesky_test() failed.");
  
  result = convertSpectralCholesky_test();
  if (!result) error("convertSpectralCholesky_test() failed.");
  
  result = convertSTToCholesky_test();
  if (!result) error("convertSTToCholesky_test() failed.");
  
  result = convertPriorCorrelationToCholesky_test();
  if (!result) error("convertPriorCorrelationToCholesky_test() failed.");
  
  result = convertCovarianceToPriorCorrelation_test();
  if (!result) error("convertCovarianceToPriorCorrelation_test() failed.");
  
  result = convertSDCorrelationToPriorCorrelation_test();
  if (!result) error("convertSDCorrelationToPriorCorrelation_test() failed.");
  
  result = convertSpectralToPriorCorrelation_test();
  if (!result) error("convertSpectralToPriorCorrelation_test() failed.");
  
  result = convertSTToPriorCorrelation_test();
  if (!result) error("convertSTToPriorCorrelation_test() failed.");
  
  result = convertCholeskyToPriorCorrelation_test();
  if (!result) error("convertCholeskyToPriorCorrelation_test() failed.");
  
  return (ScalarLogical(TRUE));
}

static int convertSDCorrelationToCovariance_test()
{
  int arrayLength = covarianceDim * covarianceDim;
  double result[arrayLength];
  convertSDCorrelationToCovariance(sdCorrelationArray, sdCorrelationDim, result);
  
  return(allApproximatelyEqual(result, covarianceMatrix, arrayLength, TEST_TOLERANCE));
}

static int convertSpectralToCovariance_test()
{
  int arrayLength = covarianceDim * covarianceDim;
  double result[arrayLength];
  convertSpectralToCovariance(spectralArray, spectralDim, result);
  
  return(allApproximatelyEqual(result, covarianceMatrix, arrayLength, TEST_TOLERANCE));
}

static int convertSTToCovariance_test()
{
  int arrayLength = covarianceDim * covarianceDim;
  double result[arrayLength];
  convertSTToCovariance(stArray, stDim, result);
  
  return(allApproximatelyEqual(result, covarianceMatrix, arrayLength, TEST_TOLERANCE));
}

static int convertSTToCovarianceWithScale_test()
{
  const double correctAnswer[] = {    4.1713644054497, -0.290425688231175, -0.065876554764186,
                                   -0.290425688231175,  0.598077272138685,  -0.56634425604058,
                                   -0.065876554764186,  -0.56634425604058,   3.93976032987112  };

  int arrayLength = covarianceDim * covarianceDim;
  double result[arrayLength];
  convertSTToCovarianceWithScale(stArray, stDim, result, 0.5);
  
  return(allApproximatelyEqual(result, correctAnswer, arrayLength, TEST_TOLERANCE));
}

static int convertCholeskyToCovariance_test()
{
  int arrayLength = covarianceDim * covarianceDim;
  double result[arrayLength];
  convertCholeskyToCovariance(choleskyMatrix, choleskyDim, result);
  
  return(allApproximatelyEqual(result, covarianceMatrix, arrayLength, TEST_TOLERANCE));
}

static int convertPriorCorrelationToCovariance_test()
{
  int arrayLength = covarianceDim * covarianceDim;
  double result[arrayLength];
  convertPriorCorrelationToCovariance(priorCorrelationArray, priorCorrelationDim, result);
  
  return(allApproximatelyEqual(result, covarianceMatrix, arrayLength, TEST_TOLERANCE));
}

static int convertSTToCovarianceInverse_test()
{
  int arrayLength = covarianceInverseDim * covarianceInverseDim;
  double result[arrayLength];
  convertSTToCovarianceInverse(stArray, stDim, result);
  
  return(allApproximatelyEqual(result, covarianceMatrixInverse,
                               arrayLength, TEST_TOLERANCE));
}

static int convertSTToCovarianceInverseWithScale_test()
{
  const double correctAnswer[] = { 0.0625595876778298, 0.0363124038521897, 0.0062659982794447,
                                   0.0363124038521897,  0.504950116245001, 0.0731942326026723,
                                   0.0062659982794447, 0.0731942326026723, 0.0740821499662589  };
                                   
  int arrayLength = covarianceInverseDim * covarianceInverseDim;
  double result[arrayLength];
  convertSTToCovarianceInverseWithScale(stArray, stDim, result, 0.5);
  
  return(allApproximatelyEqual(result, correctAnswer, arrayLength, TEST_TOLERANCE));
}



static int convertCovarianceToSDCorrelation_test() {
  int arrayLength = sdCorrelationDim + sdCorrelationDim * (sdCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertCovarianceToSDCorrelation(covarianceMatrix, covarianceDim, result);
  
  return(allApproximatelyEqual(result, sdCorrelationArray, arrayLength, TEST_TOLERANCE));
}

static int convertSpectralToSDCorrelation_test() {
  int arrayLength = sdCorrelationDim + sdCorrelationDim * (sdCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertSpectralToSDCorrelation(spectralArray, spectralDim, result);
  
  return(allApproximatelyEqual(result, sdCorrelationArray, arrayLength, TEST_TOLERANCE));
}

static int convertSTToSDCorrelation_test() {
  int arrayLength = sdCorrelationDim + sdCorrelationDim * (sdCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertSTToSDCorrelation(stArray, stDim, result);
  
  return(allApproximatelyEqual(result, sdCorrelationArray, arrayLength, TEST_TOLERANCE));
}

static int convertCholeskyToSDCorrelation_test() {
  int arrayLength = sdCorrelationDim + sdCorrelationDim * (sdCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertCholeskyToSDCorrelation(choleskyMatrix, choleskyDim, result);
  
  return(allApproximatelyEqual(result, sdCorrelationArray, arrayLength, TEST_TOLERANCE));
}

static int convertPriorCorrelationToSDCorrelation_test() {
  int arrayLength = sdCorrelationDim + sdCorrelationDim * (sdCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertPriorCorrelationToSDCorrelation(priorCorrelationArray, priorCorrelationDim, result);
  
  return(allApproximatelyEqual(result, sdCorrelationArray, arrayLength, TEST_TOLERANCE));
}


static int convertCovarianceToSpectral_test() {
  int arrayLength = covarianceDim * covarianceDim;
  double testMatrix[arrayLength];
  for (int i = 0; i < arrayLength; ++i) testMatrix[i] = covarianceMatrix[i];
  
  arrayLength = spectralDim + spectralDim * (spectralDim - 1) / 2;
  double result[arrayLength];
  convertCovarianceToSpectral(testMatrix, covarianceDim, result);
  
  return(allApproximatelyEqual(result, spectralArray, arrayLength, TEST_TOLERANCE));
}

static int convertSDCorrelationToSpectral_test() {
  int arrayLength = spectralDim + spectralDim * (spectralDim - 1) / 2;
  double result[arrayLength];
  convertSDCorrelationToSpectral(sdCorrelationArray, sdCorrelationDim, result);
  
  return(allApproximatelyEqual(result, spectralArray, spectralDim, TEST_TOLERANCE));
}

static int convertSTToSpectral_test() {
  int arrayLength = spectralDim + spectralDim * (spectralDim - 1) / 2;
  double result[arrayLength];
  convertSTToSpectral(stArray, stDim, result);
  
  return(allApproximatelyEqual(result, spectralArray, arrayLength, TEST_TOLERANCE));
}

static int convertCholeskyToSpectral_test() {
  int arrayLength = spectralDim + spectralDim * (spectralDim - 1) / 2;
  double result[arrayLength];
  convertCholeskyToSpectral(choleskyMatrix, choleskyDim, result);
  
  return(allApproximatelyEqual(result, spectralArray, arrayLength, TEST_TOLERANCE));
}

static int convertPriorCorrelationToSpectral_test() {
  int arrayLength = spectralDim + spectralDim * (spectralDim - 1) / 2;
  double result[arrayLength];
  convertPriorCorrelationToSpectral(priorCorrelationArray, priorCorrelationDim, result);
  
  return(allApproximatelyEqual(result, spectralArray, arrayLength, TEST_TOLERANCE));
}



static int convertCovarianceToST_test() {
  int arrayLength = covarianceDim * covarianceDim;
  double testMatrix[arrayLength];
  for (int i = 0; i < arrayLength; ++i) testMatrix[i] = covarianceMatrix[i];
  
  arrayLength = stDim + stDim * (stDim - 1) / 2;
  double result[arrayLength];
  convertCovarianceToST(testMatrix, covarianceDim, result);
  
  return(allApproximatelyEqual(result, stArray, arrayLength, TEST_TOLERANCE));
}

static int convertSDCorrelationToST_test() {
  int arrayLength = stDim + stDim * (stDim - 1) / 2;
  double result[arrayLength];
  convertSDCorrelationToST(sdCorrelationArray, sdCorrelationDim, result);
    
  return(allApproximatelyEqual(result, stArray, arrayLength, TEST_TOLERANCE));
}

static int convertSpectralToST_test() {
  int arrayLength = stDim + stDim * (stDim - 1) / 2;
  double result[arrayLength];
  convertSpectralToST(spectralArray, spectralDim, result);
  
  return(allApproximatelyEqual(result, stArray, arrayLength, TEST_TOLERANCE));
}

static int convertCholeskyToST_test() {
  int arrayLength = stDim + stDim * (stDim - 1) / 2;
  double result[arrayLength];
  convertCholeskyToST(choleskyMatrix, choleskyDim, result);
  
  return(allApproximatelyEqual(result, stArray, arrayLength, TEST_TOLERANCE));
}

static int convertPriorCorrelationToST_test() {
  int arrayLength = stDim + stDim * (stDim - 1) / 2;
  double result[arrayLength];
  convertPriorCorrelationToST(priorCorrelationArray, priorCorrelationDim, result);
  
  return(allApproximatelyEqual(result, stArray, arrayLength, TEST_TOLERANCE));
}


static int convertCovarianceToCholesky_test() {
  int arrayLength = choleskyDim * choleskyDim;
  double result[arrayLength];
  convertCovarianceToCholesky(covarianceMatrix, covarianceDim, result);
  
  return(allApproximatelyEqual(result, choleskyMatrix, arrayLength, TEST_TOLERANCE));
}
static int convertSDCorrelationToCholesky_test() {
  int arrayLength = choleskyDim * choleskyDim;
  double result[arrayLength];
  convertSDCorrelationToCholesky(sdCorrelationArray, sdCorrelationDim, result);
  
  return(allApproximatelyEqual(result, choleskyMatrix, arrayLength, TEST_TOLERANCE));
}

static int convertSpectralCholesky_test() {
  int arrayLength = choleskyDim * choleskyDim;
  double result[arrayLength];
  convertSpectralToCholesky(spectralArray, spectralDim, result);
  
  return(allApproximatelyEqual(result, choleskyMatrix, arrayLength, TEST_TOLERANCE));
}

static int convertSTToCholesky_test() {
  int arrayLength = choleskyDim * choleskyDim;
  double result[arrayLength];
  convertSTToCholesky(stArray, stDim, result);
  
  return(allApproximatelyEqual(result, choleskyMatrix, arrayLength, TEST_TOLERANCE));
}

static int convertPriorCorrelationToCholesky_test() {
  int arrayLength = choleskyDim * choleskyDim;
  double result[arrayLength];
  convertPriorCorrelationToCholesky(priorCorrelationArray, priorCorrelationDim, result);
  
  return(allApproximatelyEqual(result, choleskyMatrix, arrayLength, TEST_TOLERANCE));
}

static int convertCovarianceToPriorCorrelation_test() {
  double testMatrix[covarianceDim * covarianceDim];
  for (int i = 0; i < covarianceDim * covarianceDim; ++i) testMatrix[i] = covarianceMatrix[i];
  
  int arrayLength = priorCorrelationDim + priorCorrelationDim + priorCorrelationDim * (priorCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertCovarianceToPriorCorrelation(testMatrix, covarianceDim, result);
  
  return(allApproximatelyEqual(result, priorCorrelationArray, arrayLength, TEST_TOLERANCE));
}

static int convertSDCorrelationToPriorCorrelation_test() {
  int arrayLength = priorCorrelationDim + priorCorrelationDim + priorCorrelationDim * (priorCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertSDCorrelationToPriorCorrelation(sdCorrelationArray, sdCorrelationDim, result);
  
  return(allApproximatelyEqual(result, priorCorrelationArray, arrayLength, TEST_TOLERANCE));
}

static int convertSpectralToPriorCorrelation_test() {
  int arrayLength = priorCorrelationDim + priorCorrelationDim + priorCorrelationDim * (priorCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertSpectralToPriorCorrelation(spectralArray, spectralDim, result);
  
  return(allApproximatelyEqual(result, priorCorrelationArray, arrayLength, TEST_TOLERANCE));
}

static int convertSTToPriorCorrelation_test() {
  int arrayLength = priorCorrelationDim + priorCorrelationDim + priorCorrelationDim * (priorCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertSTToPriorCorrelation(stArray, stDim, result);
  
  return(allApproximatelyEqual(result, priorCorrelationArray, arrayLength, TEST_TOLERANCE));
}

static int convertCholeskyToPriorCorrelation_test() {
  int arrayLength = priorCorrelationDim + priorCorrelationDim + priorCorrelationDim * (priorCorrelationDim - 1) / 2;
  double result[arrayLength];
  convertCholeskyToPriorCorrelation(choleskyMatrix, choleskyDim, result);
  
  return(allApproximatelyEqual(result, priorCorrelationArray, arrayLength, TEST_TOLERANCE));
}
