// This file should contain utility functions that are used in
// lmm, glmm, and nlmm models.

#ifndef BLME_LMER_COMMON_H
#define BLME_LMER_COMMON_H

#include <R.h>
#include <Rinternals.h>
#include <Matrix.h> // cholmod definitions

typedef struct {
  int *factorDimensions;
  int *numGroupsPerFactor;
  int maxFactorDimension;
} SparseMatrixStructure;

void
getSparseContentAndStructure(const SEXP stExpression, const int *sparseRowForFactor,
                             double **stMatrices, SparseMatrixStructure *sparseMatrixStructure);

/**
 * Return the index of the factor associated with parameter index ind
 *
 * @param ind an index in [0, Gp[nt] - 1]
 * @param nt total number of terms
 * @param Gp group pointers, a vector of length nt+1 with Gp[0] = 0
 *
 * @return idnex of assicated factor
 */
int getFactorForSparseRow(int row, int numFactors,
                          const int *sparseRowsForFactor);

/**
 * Permute the vector src according to perm into dest
 *
 * @param dest destination
 * @param src source
 * @param perm NULL or 0-based permutation of length n
 * @param n length of src, dest and perm
 *
 * @return dest
 *
 * \note If perm is NULL the first n elements of src are copied to dest.
 */
double*
applyPermutation(double *destination, const double *source, const int *permutation, int numValues);

// dest = PAX, A is sparse, X is dense and the
// permutation is specified by an integer array
void
multiplyWithPermutation(double *destination, const int *permutation,
                        const CHM_SP sparseDesignMatrix,
                        const double *denseDesignMatrix, int numDenseColumns);

/**
 * Update the sparse model matrix A from the Zt and ST slots, where
 * A = S'T'Z'.
 *
 * @param regression an mer object
 */
void rotateSparseDesignMatrix(SEXP regression);

#endif // BLME_LMER_COMMON_H
