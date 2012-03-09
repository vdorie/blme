#include <R_ext/Lapack.h>        /* for Lapack and BLAS */

#include <R.h>
#include <Rdefines.h>

#include "matrix.h"
#include "util.h"

extern cholmod_common cholmodCommon;

void multiplyMatrices(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z)
{
  char *transa = "N", *transb = "N";
  int i,  j, k;
  double one = 1.0, zero = 0.0;
  long double sum;
  Rboolean have_na = FALSE;
  
  if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
    /* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
     * The test is only O(n) here
     */
    for (i = 0; i < nrx*ncx; i++)
	    if (ISNAN(x[i])) {have_na = TRUE; break;}
    if (!have_na)
	    for (i = 0; i < nry*ncy; i++)
        if (ISNAN(y[i])) {have_na = TRUE; break;}
    if (have_na) {
	    for (i = 0; i < nrx; i++)
        for (k = 0; k < ncy; k++) {
          sum = 0.0;
          for (j = 0; j < ncx; j++)
            sum += x[i + j * nrx] * y[j + k * nry];
          z[i + k * nrx] = sum;
        }
    } else
	    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
                      x, &nrx, y, &nry, &zero, z, &nrx);
  } else /* zero-extent operations should return zeroes */
    for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}

void multiplyMatricesWithTranspose(double *x, int nrx, int ncx, int transposeX,
                                   double *y, int nry, int ncy, int transposeY,
                                   double *z)
{
  int i, j, k;
  double one = 1.0, zero = 0.0;
  long double sum;
  Rboolean have_na = FALSE;
  
  if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
    /* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
     * The test is only O(n) here
     */
    for (i = 0; i < nrx*ncx; i++)
	    if (ISNAN(x[i])) {have_na = TRUE; break;}
    if (!have_na)
	    for (i = 0; i < nry*ncy; i++)
        if (ISNAN(y[i])) {have_na = TRUE; break;}
    if (have_na) {
	    for (i = 0; i < nrx; i++)
        for (k = 0; k < ncy; k++) {
          sum = 0.0;
          for (j = 0; j < ncx; j++)
            sum += x[i + j * nrx] * y[j + k * nry];
          z[i + k * nrx] = sum;
        }
    } else {
      char *transX = (transposeX ? "T" : "N");
      char *transY = (transposeY ? "T" : "N");
      int m = (transposeX ? ncx : nrx);
      int n = (transposeY ? nry : ncy);
          k = (transposeX ? nrx : ncx);
      
	    F77_CALL(dgemm)(transX, transY, &m, &n, &k,
                      &one, x, &nrx, y, &nry, &zero, z, &m);
    }
  } else /* zero-extent operations should return zeroes */
    for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}

void singleMatrixCrossproduct(const double *source, int numRows, int numCols,
                              double *target, int useTranspose, triangleType_t triangleType)
{
  // for us, transpose is AA', for BLAS it is the reverse (A'A)
  char *shouldTransposeMatrix = (useTranspose ? "N" : "T");
  char *useLowerTriangle = (triangleType == TRIANGLE_TYPE_UPPER ? "U" : "L");
  
  int outputDimension  = (useTranspose ? numRows : numCols);
  int inputDimension   = (useTranspose ? numCols : numRows); // if user specifies transpose, he/she wants AA'
  double d_one[]  = { 1.0, 0 };
  double d_zero[] = { 0.0, 0 };
  
  F77_CALL(dsyrk)(useLowerTriangle, shouldTransposeMatrix,
                  &outputDimension,
                  &inputDimension, d_one, source,
                  &numRows, d_zero,
                  target, &outputDimension);
  
  if (triangleType != TRIANGLE_TYPE_BOTH) return;
  
  // copy in rest of product
  for (int col = 1; col < outputDimension; ++col) {
    for (int row = 0; row < col; ++row) target[row + outputDimension * col] = target[col + outputDimension * row];
  }
}

void singleTriangularMatrixCrossproduct(const double *source, int dim, int useTranspose,
                                        triangleType_t sourceTriangleType, double *target)
{
  int offset;
  double value, diagonalElement;
  if (useTranspose == TRUE && sourceTriangleType == TRIANGLE_TYPE_LOWER) {
    for (int col = 0; col < dim; ++col) {
      offset = col * (dim + 1);
      
      // compute diagonal
      diagonalElement = source[offset];
      
      target[offset] = diagonalElement * diagonalElement;
      for (int colDecrement = 1; colDecrement <= col; ++colDecrement) {
        value = source[offset - colDecrement * dim];
        
        target[offset] += value * value;
      }
      ++offset;
      
      // compute lower triangle
      for (int row = col + 1; row < dim; ++row) {
        target[offset] = diagonalElement * source[offset];
        
        for (int k = 0; k <= col - 1; ++k) {
          target[offset] += source[row + k * dim] * source[col + k * dim];
        }
        
        
        // copy into upper triangle
        target[col + row * dim] = target[offset];
        
        ++offset;
      }
    }
  } else if (useTranspose == FALSE && sourceTriangleType == TRIANGLE_TYPE_UPPER) {
    for (int row = 0; row < dim; ++row) {
      offset = row * (dim + 1);
      
      // compute diagonal
      diagonalElement = source[offset];
      
      target[offset] = diagonalElement * diagonalElement;
      for (int rowDecrement = 1; rowDecrement <= row; ++rowDecrement) {
        value = source[offset - rowDecrement];
        
        target[offset] += value * value;
      }
      offset += dim;
      
      // compute upper triangle
      for (int col = row + 1; col < dim; ++col) {
        target[offset] = diagonalElement * source[offset];
        
        for (int k = 0; k < row - 1; ++k) {
          target[offset] += source[row + k * dim] * source[col + k * dim];
        }
        
        
        // copy into lower triangle
        target[row + col * dim] = target[offset];
        
        offset += dim;
      }
    }
  } else if (useTranspose == FALSE && sourceTriangleType == TRIANGLE_TYPE_LOWER) {
    Memcpy(target, (double* const) source, dim * dim);
    double d_one = 1.0;
    F77_CALL(dtrmm)("L", "L", "T", "N", &dim, &dim, &d_one, source,
                    &dim, target, &dim);

    for (int col = 1; col < dim; ++col) {
      for (int row = 0; row < col; ++row) target[row + dim * col] = target[col + dim * row];
    }
  } else {
    // UU'
    Memcpy(target, (double* const) source, dim * dim);
    double d_one = 1.0;
    F77_CALL(dtrmm)("R", "U", "T", "N", &dim, &dim, &d_one, source,
                    &dim, target, &dim);

    for (int row = 1; row < dim; ++row) {
      for (int col = 0; col < row; ++col) target[col + dim * row] = target[row + dim * col];
    }
  }
}

void singleMatrixCrossproductWithUpdate(const double *source, int numRows, int numCols,
                                        double crossproductScale, double *target,
                                        int useTranspose, triangleType_t triangleType)
{
  char *shouldTransposeMatrix = (useTranspose ? "N" : "T");
  char *useLowerTriangle = (triangleType == TRIANGLE_TYPE_UPPER ? "U" : "L");
  
  int outputDimension  = (useTranspose ? numRows : numCols);
  int inputDimension   = (useTranspose ? numCols : numRows); // if user specifies transpose, he/she wants AA'
  double d_one[]  = { 1.0, 0 };
  
  F77_CALL(dsyrk)(useLowerTriangle, shouldTransposeMatrix,
                  &outputDimension,
                  &inputDimension, &crossproductScale, source,
                  &numRows, d_one,
                  target, &outputDimension);
                  
  if (triangleType != TRIANGLE_TYPE_BOTH) return;
  
  // copy in rest of product
  for (int col = 1; col < outputDimension; ++col) {
    for (int row = 0; row < col; ++row) target[row + outputDimension * col] = target[col + outputDimension * row];
  }
}

void applyMatrixToVector(const double *sourceMatrix, int numRows, int numCols, int useTranspose,
                            const double *sourceVector, double *target)
{
  char *shouldTransposeMatrix = (useTranspose ? "T" : "N");
  
  int i_one = 1;
  double d_one  = 1.0;
  double d_zero = 0.0;
  
  F77_CALL(dgemv)(shouldTransposeMatrix,
                  &numRows, &numCols, &d_one,
                  sourceMatrix, &numRows,
                  sourceVector,
                  &i_one, &d_zero,
                  target,
                  &i_one);
}
void applyMatrixToVectorWithUpdate(const double *sourceMatrix, int numRows, int numCols, int useTranspose,
                                      const double *sourceVector, double productScale, double *target)
{
  char *shouldTransposeMatrix = (useTranspose ? "T" : "N");
  
  int i_one = 1;
  double d_one = 1.0;
  
  F77_CALL(dgemv)(shouldTransposeMatrix,
                  &numRows, &numCols, &productScale,
                  sourceMatrix,
                  &numRows,
                  sourceVector,
                  &i_one, &d_one,
                  target,
                  &i_one);
}

void transposeMatrix(const double *source, int sourceNumRows, int sourceNumCols,
                     double *target)
{
  int sourceOffset = 0;
  int targetOffset;
  
  for (int sourceCol = 0; sourceCol < sourceNumCols; ++sourceCol) {
    targetOffset = sourceCol;
    for (int sourceRow = 0; sourceRow < sourceNumRows; ++sourceRow) {
      target[targetOffset] = source[sourceOffset++];
      targetOffset += sourceNumCols;
    }
  }
}

// eigen values and eigen vectors from a positive definite, symmetric matrix
// source is dim x dim, values dim x 1, vectors dim x dim
// clobbers the source matrix;
int getSpectralDecompositionOfPositiveDefiniteMatrix(double *source, int dim,
                                                      double *values, double *vectors)
{
  double tolerance = 0.0;
  
  int numEigenvaluesFound;
  
  double optimalDoubleWorkLength = 0.0;
  int doubleWorkLength = -1;
  int optimalIntegerWorkLength = 0;
  int integerWorkLength = -1;
  
  int integerSupport[2 * dim];
  int lapackResult;
  
  double vl = 0.0, vu = 0.0; // not useful
  int    il = 0,   iu = 0;
  
  F77_CALL(dsyevr)("V", "A", "L",
                   &dim, source, &dim,
                   &vl, &vu, &il, &iu,
                   &tolerance, &numEigenvaluesFound,
                   values, vectors, &dim,
                   integerSupport,
                   &optimalDoubleWorkLength, &doubleWorkLength,
                   &optimalIntegerWorkLength, &integerWorkLength,
                   &lapackResult);
    
  if (lapackResult != 0) return(lapackResult);
  
  doubleWorkLength = (int) (optimalDoubleWorkLength + 0.5);
  double *doubleWork = (double *) alloca(doubleWorkLength * sizeof(double));
  
  integerWorkLength = optimalIntegerWorkLength;
  int *integerWork = (int *) alloca(integerWorkLength * sizeof(int));
  R_CheckStack();
  
  F77_CALL(dsyevr)("V", "A", "L", // compute 'v'ectors, 'a'll of them, use 'l'ower triangle of source
                   &dim, source, &dim,
                   &vl, &vu, &il, &iu,
                   &tolerance, &numEigenvaluesFound,
                   values, vectors, &dim,
                   integerSupport,
                   doubleWork, &doubleWorkLength,
                   integerWork, &integerWorkLength,
                   &lapackResult);
  
  if (lapackResult != 0) return(lapackResult);

  // swap the ordering so the eigenvalues and vectors are in descending order
  double tempValue;
  double tempVector[dim];
  
  for (int i = 0; i < dim / 2; ++i) {
    Memcpy(tempVector, vectors + i * dim, dim);
    Memcpy(vectors + i * dim, vectors + (dim - i - 1) * dim, dim);
    Memcpy(vectors + (dim - i - 1) * dim, tempVector, dim);
    tempValue = values[i];
    values[i] = values[dim - i - 1];
    values[dim - i - 1] = tempValue;
  }

  return (lapackResult);
}

// specified triangle of real, symmetric, positive definite matrix
int getDenseCholeskyDecomposition(double *matrix, int dim, triangleType_t triangleType)
{
  if (triangleType == TRIANGLE_TYPE_BOTH) error("\"Both\" specified as triangle type for a cholesky decomposition.");
  
  char *useLowerTriangle = (triangleType == TRIANGLE_TYPE_UPPER ? "U" : "L");
  
  int lapackResult;
  // find the Cholesky factor of the source
  F77_CALL(dpotrf)(useLowerTriangle, &dim, matrix, &dim, &lapackResult);
  
  return(lapackResult);
}

// mostly a code rip of what is in Lapack.c. This is only here since I can't
// readily find a C only version
double getLogDeterminant(const double *matrix, int dim, int *sign)
{
  double modulus;
  
  double *decomposition = (double *) alloca(dim * dim * sizeof(double));
  int *pivotIndices = (int *) alloca(dim * sizeof(int));
  R_CheckStack();
  
  Memcpy(decomposition, matrix, dim * dim);
  
  *sign = 1;
  
  int lapackResult;
  // take the LU decomposition.
  // the diagonal of L is unit, and ignored. The diagonal of U is then
  // the diagonal of the result, and the determinant is the determined
  // by the diagonal of U up to the sign, which is given by the pivot
  F77_CALL(dgetrf)(&dim, &dim, decomposition, &dim, pivotIndices, &lapackResult);
  if (lapackResult < 0) {
    error("error taking LU decomposition, argument %d is illegal", -lapackResult);
  } else if (lapackResult > 0) {
    modulus = (ISNAN(matrix[(lapackResult - 1) * dim + lapackResult - 1]) ? R_NaN : R_NegInf);
    return (modulus);
  }
  
  for (int i = 0; i < dim; ++i) {
    if (pivotIndices[i] != (i + 1)) *sign = -*sign;
  }
  
  modulus = 0.0;
  for (int i = 0; i < dim; ++i) {
    double diagonalElement = decomposition[i * (dim + 1)];
    
    modulus += log(diagonalElement < 0 ? -diagonalElement : diagonalElement);
    if (diagonalElement < 0) *sign = -*sign;
  }
  
  return (modulus);
}

double getLogDeterminantOfTriangularMatrix(const double *matrix, int dim, int *sign)
{
  *sign = 1;
  
  double modulus = 0.0;
  for (int i = 0; i < dim; ++i) {
    double diagonalElement = matrix[i * (dim + 1)];
    
    modulus += log(diagonalElement < 0 ? -diagonalElement : diagonalElement);
    if (diagonalElement < 0) *sign = -*sign;
  }
  
  return (modulus);
}


double getLogDeterminantOfPositiveDefiniteMatrix(const double *matrix, int dim)
{
  double *decomposition = (double *) alloca(dim * dim * sizeof(double));
  R_CheckStack();
  
  Memcpy(decomposition, matrix, dim * dim);
  
  int lapackResult;
  // take the Cholesky factorization. Since LL' = A, det(A) = 2 det(L) = 2 * prod * diag(L)
  F77_CALL(dpotrf)("L", &dim, decomposition, &dim, &lapackResult);
  if (lapackResult != 0) {
    if (lapackResult < 0) {
      error("error taking cholesky factorization, argument %d is invalid", -lapackResult);
    } else {
      return (R_NegInf);
    }
  }
  
  int sign;
  double result;
  result = 2.0 * getLogDeterminantOfTriangularMatrix(decomposition, dim, &sign);
  return (result);
}

// AX = B, X = A^{-1}B; A: lhs x lhs, B: lhs x rhs; X: lhs x rhs
int solveSystem(double *leftHandSide, int commonDim,
                double *rightHandSide, int numColsRightHandSide,
                double *solution)
{
  double factor[commonDim * commonDim];
  int factorPivotIndices[commonDim];
  double rowScaleFactors[commonDim];
  double colScaleFactors[commonDim];
  
  char equilibration;

  double reciprocalConditionNumber;
  double forwardErrorBound[numColsRightHandSide];
  double backwardErrorBound[numColsRightHandSide];
  
  double doubleScratch[4 * commonDim];
  int    integerScratch[commonDim];
  
  int lapackResult;

  int factorize = (int) 'E'; // spec is a char, but R has it as an int

  F77_CALL(dgesvx)(&factorize, "N", &commonDim, &numColsRightHandSide,
                   leftHandSide, &commonDim, factor, &commonDim,
                   factorPivotIndices, &equilibration, rowScaleFactors,
                   colScaleFactors, rightHandSide, &commonDim,
                   solution, &commonDim, &reciprocalConditionNumber, forwardErrorBound,
                   backwardErrorBound, doubleScratch, integerScratch, &lapackResult);

  return (lapackResult);
}

// same as solve, but A is assumed to be symmetric
// only the lower triangle is used, solution stored in rightHandSide matrix
int solveSymmetricSystem(double *leftHandSide, int commonDim,
                         double *rightHandSide, int numColsRightHandSide)
{
  
  int lapackResult;
  F77_CALL(dposv)("L", &commonDim, &numColsRightHandSide,
                  leftHandSide, &commonDim, rightHandSide,
                  &commonDim, &lapackResult);
                  
  return(lapackResult);
}

void solveSparseCholeskySystem(int operationType, const_CHM_FR leftHandSide,
                               const double *rightHandSide, int numColsRightHandSide,
                               double *target)
{
  int commonDim = leftHandSide->n;
  
  CHM_DN chol_rightHandSide = N_AS_CHM_DN((double *) rightHandSide, commonDim, numColsRightHandSide);
  CHM_DN solution;
  R_CheckStack();
  
  solution = M_cholmod_solve(operationType,
                             leftHandSide,
                             chol_rightHandSide,
                             &cholmodCommon);
  if (!solution) error("cholmod_solve (%d) failed", operationType);
  
  // copy result into the dense, non-cholmod matrix
  Memcpy(target, (const double *) (solution->x),
         commonDim * numColsRightHandSide);
  M_cholmod_free_dense(&solution, &cholmodCommon);
}

void invertLowerTriangularMatrix(const double *source, int dim, double *target)
{
  double diagonal;
  
  int offset = 0;
  for (int col = 0; col < dim; ++col) {
    // 0 out the upper triangle
    for (int row = 0; row < col; ++row) target[offset++] = 0.0;
    
    diagonal = 1.0 / source[offset];
    target[offset++] = diagonal;
    
    for (int row = col + 1; row < dim; ++row) {
      target[offset] = -source[offset] * diagonal;
      
      for (int k = col + 1; k < row; ++k) {
        target[offset] -= target[k + col * dim] * source[row + k * dim];
      }
      target[offset++] /= source[row * (dim + 1)];
    }
  }
}

void invertUnitLowerTriangularMatrix(const double *source, int dim, double *target)
{
  int offset = 0;
  for (int col = 0; col < dim; ++col) {
    // 0 out the upper triangle
    for (int row = 0; row < col; ++row) target[offset++] = 0.0;
    
    target[offset++] = 1.0;
    
    for (int row = col + 1; row < dim; ++row) {
      target[offset] = -source[offset];
      
      for (int k = col + 1; k < row; ++k) {
        target[offset] -= target[k + col * dim] * source[row + k * dim];
      }
      
      ++offset;
    }
  }
}

void invertUpperTriangularMatrix(const double *source, int dim, double *target)
{
  double diagonal;
  
  int offset = dim * dim - 1;
  for (int col = dim - 1; col >= 0; --col) {
    // 0 out the lower triangle
    for (int row = dim - 1; row > col; --row) {
      target[offset--] = 0.0;
    }
    
    diagonal = 1.0 / source[offset];
    target[offset--] = diagonal;
    
    for (int row = col - 1; row >= 0; --row) {
      target[offset] = -source[offset] * diagonal;
      
      for (int k = col - 1; k > row; --k) {
        target[offset] -= target[k + col * dim] * source[row + k * dim];
      }
      target[offset--] /= source[row * (dim + 1)];
    }
  }
}

void invertUnitUpperTriangularMatrix(const double *source, int dim, double *target)
{
  int offset = dim * dim - 1;
  for (int col = dim - 1; col >= 0; --col) {
    // 0 out the lower triangle
    for (int row = dim - 1; row > col; --row) {
      target[offset--] = 0.0;
    }
    
    target[offset--] = 1.0;
    
    for (int row = col - 1; row >= 0; --row) {
      target[offset] = -source[offset];
      
      for (int k = col - 1; k > row; --k) {
        target[offset] -= target[k + col * dim] * source[row + k * dim];
      }
      --offset;
    }
  }
}

// Q = (I - A)(I + A)^{-1}
void getCayleyTransform(const double *source, int dim, double *target)
{
  double temp1[dim * dim];
  double temp2[dim * dim];
  double temp3[dim * dim];
  
  double *plusMatrix     = temp1;          // temp1 in use; I+A
  double *minusMatrix    = temp2;          // temp2 in use; I-A
  double *identity       = target;         //               I
  
  // copy in I+A to lower triangle for inverse, and set scratch to I,
  // copy I-A to temp
  int offset = 0;
  for (int col = 0; col < dim; ++col) {
    offset = col * (dim + 1);
    
    identity[offset] = 1.0;
    plusMatrix[offset] = 1.0;
    minusMatrix[offset] = 1.0;
    ++offset;

    for (int row = col + 1; row < dim; ++row) {
      identity[offset] = identity[col + row * dim] = 0.0;
      
      plusMatrix[offset] = source[offset];
      plusMatrix[col + row * dim] = source[col + row * dim];
      
      minusMatrix[offset] = -source[offset];
      minusMatrix[col + row * dim] = -source[col + row * dim];
      
      ++offset;
    }
  }
  
  double *inverseMatrix  = temp3;           // temp3 in use; (I+Q)^{-1}
  
  int lapackResult = solveSystem(plusMatrix, dim, identity, dim, inverseMatrix);
                                            // temp1 free
  
  if (lapackResult != 0) {
    if (lapackResult < 0) {
      error("error in call to LAPACK routine 'dgesvx': argument %d illegal", -lapackResult);
    } else if (lapackResult <= dim ){
      error("error in call to LAPACK routine 'dgesvx': factor U(%d) is exactly 0 so that U is singular", lapackResult);
    } else {
      error("error in call to LAPACK routine 'dgesvx': reciprocal condition estimate below machine tolerance", lapackResult);
    }
  }
  
  multiplyMatrices(minusMatrix, dim, dim,      // temp2, temp3 free
                   inverseMatrix, dim, dim,
                   target);
}

// A = (I - Q)(I + Q)^{-1}
void getInverseCayleyTransform(const double *source, int dim, double *target)
{
  double temp1[dim * dim];
  double temp2[dim * dim];
  double temp3[dim * dim];
  
  double *plusMatrix     = temp1;          // temp1 in use; I+Q
  double *minusMatrix    = temp2;          // temp2 in use; I-Q
  double *identity       = target;         //               I
  
  
  int offset = 0;
  for (int col = 0; col < dim; ++col) {
    for (int row = 0; row < col; ++row) {
      identity[offset]    = 0.0;
      plusMatrix[offset]  = source[offset];
      minusMatrix[offset] = -source[offset];
      ++offset;
    }
    
    identity[offset]    = 1.0;
    plusMatrix[offset]  = 1.0 + source[offset];
    minusMatrix[offset] = 1.0 - source[offset];
    ++offset;
    
    for (int row = col + 1; row < dim; ++row) {
      identity[offset]    = 0.0;
      plusMatrix[offset]  = source[offset];
      minusMatrix[offset] = -source[offset];
      ++offset;
    }
  }
  
  double *inverseMatrix  = temp3;           // temp3 in use; (I+Q)^{-1}
  
  int lapackResult = solveSystem(plusMatrix, dim, identity, dim, inverseMatrix);
                                            // temp1 free
  
  if (lapackResult != 0) {
    if (lapackResult < 0) {
      error("error in call to LAPACK routine 'dgesvx': argument %d illegal", -lapackResult);
    } else if (lapackResult <= dim) {
      error("error in call to LAPACK routine 'dgesvx': factor U(%d) is exactly 0 so that U is singular", lapackResult);
    } else {
      error("error in call to LAPACK routine 'dgesvx': reciprocal condition estimate below machine tolerance", lapackResult);
    }
  }
  
  multiplyMatrices(minusMatrix, dim, dim,   // temp2, temp3 free
                   inverseMatrix, dim, dim,
                   target);
}
