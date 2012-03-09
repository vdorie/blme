#ifndef BLME_MATRIX_H
#define BLME_MATRIX_H

#include <Matrix.h> // for cholmod

// this is the matprod function from src/main/array.c, which is unfortunately not exported
// although x and y are not const, they should not change (subject to the whims of the
// Fortran implementation)
void multiplyMatrices(double *x, int numRowsX, int numColsX,
                      double *y, int numRowsY, int numColsY,
                      double *z);
                      
void multiplyMatricesWithTranspose(double *x, int numRowsX, int numColsX, int transposeX,
                                   double *y, int numRowsY, int numColsY, int transposeY,
                                   double *z);

typedef enum {
  TRIANGLE_TYPE_UPPER,
  TRIANGLE_TYPE_LOWER,
  TRIANGLE_TYPE_BOTH,
} triangleType_t;
// X'X or XX', Y is numCols x numCols or numRows x numRows depending on tranpose arg
// this is "crossprod" in R, when used on a single target.
// Note that we use R's default which is X'X instead of BLAS' view which prioritizes XX'
void singleMatrixCrossproduct(const double *source, int numRows, int numCols,
                              double *target, int useTranspose, triangleType_t triangleType);
// source' * source (default) or source * source'
void singleTriangularMatrixCrossproduct(const double *source, int dim, int useTranspose,
                                        triangleType_t sourceTriangleType, double *target);
void singleMatrixCrossproductWithUpdate(const double *source, int numRows, int numCols,
                                        double crossproductScale, double *target,
                                        int useTranspose, triangleType_t triangleType);
// y = Ax or y = A'x
void applyMatrixToVector(const double *sourceMatrix, int numRows, int numCols, int useTranpose,
                         const double *sourceVector, double *target);
// y: = Ax + b * y or y:= A'x + b*y;
void applyMatrixToVectorWithUpdate(const double *sourceMatrix, int numRows, int numCols, int useTranpose,
                                   const double *sourceVector, double productScale, double *target);
                           
void transposeMatrix(const double *source, int numRowsSource, int numColsSource,
                     double *target);


// AX = B, X = A^{-1}B; A: lhs x lhs, B: lhs x rhs; X: lhs x rhs
// possibly clobbers both input matrices
int solveSystem(double *leftHandSide, int commonDim,
                double *rightHandSide, int numColsRightHandSide,
                double *solution);
          
// same as solve, but A is assumed to be symmetric
// only the lower triangle is used, solution stored in rightHandSide matrix
int solveSymmetricSystem(double *leftHandSide, int commonDim,
                         double *rightHandSide, int numColsRightHandSide);


// see Matrix/include/cholmod.h for operation types
// wraps the construction of a dense cholmod matrix around the contents
// of RHS and copies the result from the cholmod struct to the dense
// target.
//
// If the LHS is n x n, RHS and target are n x p. Cholmod struct
// includes its dimension
void solveSparseCholeskySystem(int operationType, const_CHM_FR leftHandSide,
                               const double *rightHandSide, int numColsRightHandSide,
                               double *target);

// stores the full result, 0ing out the upper triangle
// source cannot equal target
void invertLowerTriangularMatrix(const double *source, int dim, double *target);
// stores the full result, 0ing out the upper triangle and setting
// diagonal to unity
void invertUnitLowerTriangularMatrix(const double *source, int dim,
                               double *target);
void invertUpperTriangularMatrix(const double *source, int dim, double *target);
void invertUnitUpperTriangularMatrix(const double *source, int dim,
                                     double *target);
                               

// eigen values and eigen vectors from a positive definite, symmetric matrix
// clobbers the source matrix
int getSpectralDecompositionOfPositiveDefiniteMatrix(double *source, int dim,
                                                      double *values, double *vectors);

// specified triangle of real, symmetric, positive definite matrix. leaves off
// triangle unchanged, so clean it if need-be
// result is a LAPACK error code
int getDenseCholeskyDecomposition(double *matrix, int dim, triangleType_t triangleType);

// returns the modulus. stores the sign in the given pointer
double getLogDeterminant(const double *matrix, int dim, int *sign);

double getLogDeterminantOfTriangularMatrix(const double *matrix, int dim, int *sign);
// uses a cholesky decomposition. only returns the modulus (sign must be positive).
double getLogDeterminantOfPositiveDefiniteMatrix(const double *matrix, int dim);

// Q = (I - A)(I + A)^{-1}
// A is skew-symmetric, so that only the lower triangle is necessary
void getCayleyTransform(const double *source, int dim, double *target);

// A = (I - Q)(I + Q)^{-1}
// A is skew symmetric, but the entire matrix is returned
void getInverseCayleyTransform(const double *source, int dim, double *target);

#endif // BLME_MATRIX_H
