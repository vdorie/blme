#include "parameters.h"

#include <R.h>     // for error and Memcpy
#include <Rmath.h> // for sqrt and bools

#include "matrix.h"
#include "util.h"



void convertSDCorrelationToCovariance(const double *source, int dim, double *target)
{
  if (dim == 1) {
    *target = source[0] * source[0];
    return;
  }
  
  const double *sds = source;
  const double *correlations = source + dim;
  
  for (int col = 0; col < dim; ++col) {
    target[col * (dim + 1)] = sds[col] * sds[col];
    
    for (int row = col + 1; row < dim; ++row) {
      target[row + col * dim] = target[col + row * dim] = *correlations++ * sds[col] * sds[row];
    }
  }
}

void convertSpectralToCovariance(const double *source, int dim, double *target)
{
  if (dim == 1) {
    *target = *source;
    return;
  }
  // we store the spectral decomposition, Q Lambda Q' as:
  // [ l_1 ... l_q theta_21 ... a_21 ... a_q1 ... a_q(q-1) ],
  // where a_jk are the lower triangular elements of a skew symmetric matrix
  // that are transformed into an orthogonal matrix
  
  double temp1[dim * dim];
  double temp2[dim * dim];
  
  double eigenvalue;
  const double *eigenvalues = source;
  const double *cayleyParams = source + dim;
  
  double *skewSymmetricMatrix = temp1;                  // temp1 in use
  
  int offset = 0;
  for (int col = 0; col < dim; ++col) {
    skewSymmetricMatrix[col * (dim + 1)] = 0.0;
    for (int row = col + 1; row < dim; ++row) {
      skewSymmetricMatrix[row + col * dim] = cayleyParams[offset];
      skewSymmetricMatrix[col + row * dim] = -cayleyParams[offset];
      ++offset;
    }
  }
  
  double *orthogonalMatrix = temp2;                     // temp2 in use
  getCayleyTransform(skewSymmetricMatrix, dim, orthogonalMatrix);
                                                        // temp1 free
    
  // scale rotation matrix by eigen values
  double *partialProduct = temp1;                       // temp1 in use
  offset = 0;
  for (int col = 0; col < dim; ++col) {
    eigenvalue = eigenvalues[col];
    for (int row = 0; row < dim; ++row) {
      partialProduct[offset] = orthogonalMatrix[offset] * eigenvalue;
      ++offset;
    }
  }
  
  multiplyMatricesWithTranspose(partialProduct, dim, dim, 0,
                                orthogonalMatrix, dim, dim, 1,
                                target);
}

void convertSTToCovariance(const double *source, int dim, double *target)
{
  // "ST" is lmer style, so Cov = TSST'
  // source is stored as [ s_1 ... s_q l_21 ... l_q1 ... l_q(q-1) ]; length = d + d * (d - 1) / 2, d * (d + 1) / 2
  
  const double *scales = source;
  const double *factors = source + dim;
  
  int arraySize = dim * dim;
  
  double lowerTriangle[arraySize];

  // directly scale the lower triangle so that it equals T * S
  
  double scale;
  int offset = 0;
  for (int col = 0; col < dim; ++col) {
    scale = scales[col];
    
    for (int row = 0; row < col; ++row) lowerTriangle[offset++] = 0.0;
    
    lowerTriangle[offset++] = scale;
    
    for (int row = col + 1; row < dim; ++row) {
      lowerTriangle[offset++] = scale * *factors++;
    }
  }
  
//  singleMatrixCrossproduct(lowerTriangle, dim, dim, target, TRUE, TRIANGLE_TYPE_BOTH);
  singleTriangularMatrixCrossproduct(lowerTriangle, dim, TRUE, TRIANGLE_TYPE_LOWER, target);
}

void convertSTToCovarianceWithScale(const double *source, int dim, double *target, double matrixScale)
{
  matrixScale = sqrt(matrixScale); // copy in only sqrt, since we take crossprod later
  
  // sadly is a copy/paste of above
  
  const double *scales = source;
  const double *factors = source + dim;
  
  int arraySize = dim * dim;
  
  double lowerTriangle[arraySize];

  
  double scale;
  int offset = 0;
  for (int col = 0; col < dim; ++col) {
    scale = scales[col];
    
    for (int row = 0; row < col; ++row) lowerTriangle[offset++] = 0.0;
    
    lowerTriangle[offset++] = scale * matrixScale;
    
    for (int row = col + 1; row < dim; ++row) {
      lowerTriangle[offset++] = scale * *factors++ * matrixScale;
    }
  }
  
//  singleMatrixCrossproduct(lowerTriangle, dim, dim, target, TRUE, TRIANGLE_TYPE_BOTH);
  singleTriangularMatrixCrossproduct(lowerTriangle, dim, TRUE, TRIANGLE_TYPE_LOWER, target);
}

void convertCholeskyToCovariance(const double *source, int dim, double *target)
{
//  singleMatrixCrossproduct(source, dim, dim, target, TRUE, TRIANGLE_TYPE_BOTH);
  singleTriangularMatrixCrossproduct(source, dim, TRUE, TRIANGLE_TYPE_LOWER, target);
}


void convertPriorCorrelationToCovariance(const double *source, int dim, double *target)
{
  if (dim == 1) {
    *target = source[0] * source[0] * source[1] * source[1];
    return;
  }
  
  const double *scales = source;
  
  convertSTToCovariance(source + dim, dim, target);
  
  for (int col = 0; col < dim; ++col) {
    for (int row = 0; row < dim; ++row) {
      *target++ *= scales[col] * scales[row];
    }
  }
}

void convertSTToCovarianceInverse(const double *source, int dim, double *target)
{
  // ST is lmer style, TSST'
  // source is stored as [ s_1 ... s_q l_21 ... l_q1 ... l_q(q-1) ]
  
  const double *scales = source;
  const double *factors = source + dim;
  
  int arraySize = dim * dim;
  
  double lowerTriangle[arraySize];
  double lowerTriangleInverse[arraySize];
  
  // copy into L the contents of T
  int offset = 0;
  for (int col = 0; col < dim; ++col) {
    for (int row = 0; row < col; ++row) lowerTriangle[offset++] = 0.0;
    
    lowerTriangle[offset++] = 1.0;
    
    for (int row = col + 1; row < dim; ++row) {
      lowerTriangle[offset++] = *factors++;
    }
  }
  
  invertUnitLowerTriangularMatrix((const double *) lowerTriangle, dim, lowerTriangleInverse);
  
  // TS = L, L^-1 = S^-1 T^-1
  // so we need to scale down by the inverse
  for (int row = 0; row < dim; ++row) {
    double scale = 1 / scales[row];
    
    lowerTriangleInverse[row * (dim + 1)] = scale;
    
    for (int col = 0; col < row; ++col) {
      lowerTriangleInverse[row + col * dim] *= scale;
    }
  }
  
//  singleMatrixCrossproduct(lowerTriangleInverse, dim, dim, target, FALSE, TRIANGLE_TYPE_BOTH);
  singleTriangularMatrixCrossproduct((const double *) lowerTriangleInverse, dim, FALSE, TRIANGLE_TYPE_LOWER, target);
}

void convertSTToCovarianceInverseWithScale(const double *source, int dim, double *target,
                                           double matrixScale)
{
  matrixScale = sqrt(matrixScale);
  
  // ST is lmer style, TSST'
  // source is stored as [ s_1 ... s_q l_21 ... l_q1 ... l_q(q-1) ]
  
  const double *scales = source;
  const double *factors = source + dim;
  
  int arraySize = dim * dim;
  
  double lowerTriangle[arraySize];
  double lowerTriangleInverse[arraySize];
  
  // copy into L the contents of T
  int offset = 0;
  for (int col = 0; col < dim; ++col) {
    for (int row = 0; row < col; ++row) lowerTriangle[offset++] = 0.0;
    
    lowerTriangle[offset++] = 1.0;
    
    for (int row = col + 1; row < dim; ++row) {
      lowerTriangle[offset++] = *factors++;
    }
  }
  
  invertUnitLowerTriangularMatrix((const double *) lowerTriangle, dim, lowerTriangleInverse);
  
  // TS = L, L^-1 = S^-1 T^-1
  // so we need to scale down by the inverse
  for (int row = 0; row < dim; ++row) {
    double scale = 1 / scales[row];
    
    lowerTriangleInverse[row * (dim + 1)] = scale * matrixScale;
    
    for (int col = 0; col < row; ++col) {
      lowerTriangleInverse[row + col * dim] *= scale * matrixScale;
    }
  }
  
//  singleMatrixCrossproduct(lowerTriangleInverse, dim, dim, target, FALSE, TRIANGLE_TYPE_BOTH);
  singleTriangularMatrixCrossproduct(lowerTriangleInverse, dim, FALSE, TRIANGLE_TYPE_LOWER, target);
}

void convertCovarianceToSDCorrelation(const double *source, int dim, double *target)
{
  double *targetSDs = target;
  double *targetCorrs = target + dim;
  
  int colOffset = 0;
  for (int i = 0; i < dim; ++i) {
    targetSDs[i] = sqrt(source[colOffset]);
    colOffset += dim + 1;
  }
  
  colOffset = 0;
  for (int col = 0; col < dim - 1; ++col) {
    for (int row = col + 1; row < dim; ++row) {
      *targetCorrs++ = source[colOffset + row] / (targetSDs[col] * targetSDs[row]);
    }
    colOffset += dim;
  }
}

void convertSpectralToSDCorrelation(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertSpectralToCovariance(source, dim, covarianceMatrix);
  
  convertCovarianceToSDCorrelation(covarianceMatrix, dim, target);
}

void convertSTToSDCorrelation(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertSTToCovariance(source, dim, covarianceMatrix);
    
  convertCovarianceToSDCorrelation(covarianceMatrix, dim, target);
}

void convertCholeskyToSDCorrelation(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertCholeskyToCovariance(source, dim, covarianceMatrix);
    
  convertCovarianceToSDCorrelation(covarianceMatrix, dim, target);
}

void convertPriorCorrelationToSDCorrelation(const double *source, int dim, double *target)
{
 // since the correlation prior works on XiSigmaXi, Xi diagonal, Sigma positive definite, and we
  // store this as [ xi_1 ... xi_q s_1 ... s_q t_21 .. t_q1 .. t_qq-1 ],
  // we can precompute a few offsets to make life easier
  const double *xi = source;
  const double *st = source + dim;
  
  double covarianceMatrix[dim * dim];
  
  convertSTToCovariance(st, dim, covarianceMatrix);
  
  // we are converting to a SRS, S = standard deviations, R = correlations, which we store as:
  // [ s_1 ... s_q r_21 ... r_q1 ... r_q(q-1) ]
  
  // s_i = xi_i * sqrt(sigma_ii)
  // r_ij = sigma_ij / sqrt(sigma_jj * sigma_ii) (i != j)
  
  for (int i = 0; i < dim; ++i) {
    *target++ = xi[i] * sqrt(covarianceMatrix[i * (dim + 1)]);
  }
  
  for (int col = 0; col < dim - 1; ++col) {
    double colVariance = covarianceMatrix[col * (dim + 1)];
    
    for (int row = col + 1; row < dim; ++row) {
      *target++ = covarianceMatrix[row + col * dim] / 
        sqrt(colVariance * covarianceMatrix[row * (dim + 1)]);
    }
  }
}

void convertCovarianceToSpectral(double *source, int dim, double *target)
{
  double eigenvalues[dim];
  double eigenvectors[dim * dim];
  
  int lapackResult = getSpectralDecompositionOfPositiveDefiniteMatrix(source, dim, eigenvalues, eigenvectors);
  
  if (lapackResult != 0) {
    if (lapackResult > 0) error("%d eigenvectors failed to converge, check IFAIL array.", lapackResult);
    else error("error in call to LAPACK routine 'dsyevx': argument %d illegal", -lapackResult);
  }
  
  int sign;
  getLogDeterminant(eigenvectors, dim, &sign);
  if (sign < 0) {
    // implies that the eigenvectors do not belong to SO(n), so we swap the last two columns and values
    double tempScalar = eigenvalues[dim - 2];
    eigenvalues[dim - 2] = eigenvalues[dim - 1];
    eigenvalues[dim - 1] = tempScalar;
    
    double tempVector[dim];
    Memcpy(tempVector, (const double *) eigenvectors + dim * (dim - 2), dim);
    Memcpy(eigenvectors + dim * (dim - 2), (const double *) eigenvectors + dim * (dim - 1), dim);
    Memcpy(eigenvectors + dim * (dim - 1), (const double *) tempVector, dim);
  }
    
  double skewSymmetricMatrix[dim * dim];
  
  getInverseCayleyTransform(eigenvectors, dim, skewSymmetricMatrix);
  
  for (int i = 0; i < dim; ++i) *target++ = eigenvalues[i];

  for (int col = 0; col < (dim - 1); ++col) {
    for (int row = col + 1; row < dim; ++row) *target++ = skewSymmetricMatrix[row + col * dim];
  }
}

void convertSDCorrelationToSpectral(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertSDCorrelationToCovariance(source, dim, covarianceMatrix);
  
  convertCovarianceToSpectral(covarianceMatrix, dim, target);
}

void convertSTToSpectral(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertSTToCovariance(source, dim, covarianceMatrix);
    
  convertCovarianceToSpectral(covarianceMatrix, dim, target);
}

void convertCholeskyToSpectral(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertCholeskyToCovariance(source, dim, covarianceMatrix);
  
  convertCovarianceToSpectral(covarianceMatrix, dim, target);
}

void convertPriorCorrelationToSpectral(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertPriorCorrelationToCovariance(source, dim, covarianceMatrix);
  
  convertCovarianceToSpectral(covarianceMatrix, dim, target);
}

void convertCovarianceToST(double *source, int dim, double *target)
{
  int lapackResult = getDenseCholeskyDecomposition(source, dim, TRIANGLE_TYPE_LOWER);
  
  if (lapackResult != 0) {
    if (lapackResult > 0) error("error in call to LAPACK routine 'dpotrf': leading minor of order %d is not positive definite.", lapackResult);
    else error("error in call to LAPACK routine 'dpotrf': argument %d illegal", -lapackResult);
  }
  
  convertCholeskyToST(source, dim, target);
}

void convertSDCorrelationToST(const double * source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertSDCorrelationToCovariance(source, dim, covarianceMatrix);
  
  convertCovarianceToST(covarianceMatrix, dim, target);
}

void convertSpectralToST(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertSpectralToCovariance(source, dim, covarianceMatrix);
  
  convertCovarianceToST(covarianceMatrix, dim, target);
}

void convertCholeskyToST(const double *source, int dim, double *target)
{
  // since we want Sigma=TSST', and we have L=(TS), we will have to divide
  // each column of L by its diagonal element to get T. However, we
  // ignore the fact that T has a unit diagonal and write S there instead
  
  double *scales = target;
  double *lowerTriangle = target + dim;
  
  // copy in diagonal unchanged
  int offset = 0;
  for (int i = 0; i < dim; ++i) {
    scales[i] = source[offset];
    offset += dim + 1;
  }
  
  offset = 0;
  for (int col = 0; col < dim - 1; ++col) {
    offset = col * (dim + 1) + 1;
    double scale = scales[col];
    
    for (int row = col + 1; row < dim; ++row) {
      *lowerTriangle++ = source[offset] / scale;
      ++offset;
    }
  }
}

void convertPriorCorrelationToST(const double *source, int dim, double *target)
{
  // can just push scales into diagonal of cholesky
  for (int i = 0; i < dim; ++i) {
    *target++ = source[i] * source[i + dim];
  }
  
  const double *sourceOuterScales = source;
  const double *sourceLowerTriangle = source + 2 * dim;
  
  // we have to scale each element of the lower triangle by the outer scale matrix,
  // which we essentially want to push through to the inner scale matrix.
  
  for (int col = 0; col < (dim - 1); ++col) {
    double colScale = sourceOuterScales[col];
    
    for (int row = col + 1; row < dim; ++row) {
      *target++ = *sourceLowerTriangle++ * sourceOuterScales[row] / colScale;
    }
  }
}

void convertCovarianceToCholesky(const double *source, int dim, double *target)
{
  Memcpy(target, source, dim * dim);
  
  int lapackResult = getDenseCholeskyDecomposition(target, dim, TRIANGLE_TYPE_LOWER);
  
  if (lapackResult != 0) {
    if (lapackResult > 0) error("error in call to LAPACK routine 'dpotrf': leading minor of order %d is not positive definite.", lapackResult);
    else error("error in call to LAPACK routine 'dpotrf': argument %d illegal", -lapackResult);
  }
  
  // 0 out the upper triangle
  for (int col = 1; col < dim; ++col) {
    int colOffset = col * dim;
    
    for (int row = 0; row < col; ++row) target[row + colOffset] = 0.0;
  }

}

void convertSDCorrelationToCholesky(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertSDCorrelationToCovariance(source, dim, covarianceMatrix);
  
  convertCovarianceToCholesky(covarianceMatrix, dim, target);
}

void convertSpectralToCholesky(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertSpectralToCovariance(source, dim, covarianceMatrix);
  
  convertCovarianceToCholesky(covarianceMatrix, dim, target);
}

void convertSTToCholesky(const double *source, int dim, double *target)
{
  const double *scales = source;
  const double *lowerTriangle = source + dim;
  
  for (int col = 0; col < dim; ++col) {
    for (int row = 0; row < col; ++row) *target++ = 0.0;
    
    double scale = scales[col];
    
    *target++ = scale;
    
    for (int row = col + 1; row < dim; ++row) *target++ = *lowerTriangle++ * scale;
  }
}

void convertPriorCorrelationToCholesky(const double *source, int dim, double *target)
{
  double covarianceMatrix[dim * dim];
  
  convertPriorCorrelationToCovariance(source, dim, covarianceMatrix);
  
  convertCovarianceToCholesky(covarianceMatrix, dim, target);
}


void convertCovarianceToPriorCorrelation(double *source, int dim, double *target)
{
  double st[dim * (dim + 1) / 2];
  
  convertCovarianceToST(source, dim, st);
  
  convertSTToPriorCorrelation(st, dim, target);
}

void convertSDCorrelationToPriorCorrelation(const double *source, int dim, double *target)
{
  double st[dim * (dim + 1) / 2];
  
  convertSDCorrelationToST(source, dim, st);
  
  convertSTToPriorCorrelation(st, dim, target);
}

void convertSpectralToPriorCorrelation(const double *source, int dim, double *target)
{
  double st[dim * (dim + 1) / 2];
  
  convertSpectralToST(source, dim, st);
  
  convertSTToPriorCorrelation(st, dim, target);
}

void convertSTToPriorCorrelation(const double *source, int dim, double *target)
{
  // the scale params are unidentifiable, so we split the cholesky scales by sqrt into
  // each position
  
  const double *sourceScales = source;
  const double *sourceLowerTriangle = source + dim;
  
  double *targetOuterScales = target;
  double *targetInnerScales = target + dim;
  double *targetLowerTriangle = target + 2 * dim;
  
  for (int i = 0; i < dim; ++i) {
    double scale = sqrt(sourceScales[i]);
    
    targetOuterScales[i] = scale;
    targetInnerScales[i] = scale;
  }
  
  for (int col = 0; col < (dim - 1); ++col) {
    double colScale = targetOuterScales[col];
    for (int row = col + 1; row < dim; ++row) {
      *targetLowerTriangle++ = *sourceLowerTriangle++ * colScale / targetOuterScales[row];
    }
  }
}

void convertCholeskyToPriorCorrelation(const double *source, int dim, double *target)
{
  double st[dim * (dim + 1) / 2];
  
  convertCholeskyToST(source, dim, st);
  
  convertSTToPriorCorrelation(st, dim, target);
}
