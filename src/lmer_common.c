#include "lmer_common.h"
#include "util.h"
#include "lmer.h"
#include "Syms.h"

#include <R_ext/Lapack.h>

static void
blockwiseLeftMultiplyMatrixByUpperUnitTrianglularMatrices(CHM_SP sparseDesignMatrix,
                                                          int *sparseRowForFactor, 
                                                          double **stMatrices, int numFactors,
                                                          const SparseMatrixStructure *sparseMatrixStructure);

static void
blockwiseLeftMultiplyMatrixByDiagonalMatrices(CHM_SP sparseDesignMatrix,
                                              int *sparseRowForFactor,
                                              double **stMatrices, int numFactors,
                                              const SparseMatrixStructure *sparseMatrixStructure);

/**
 * Populate the stMatrices, numColumnsPerFactor and numRowsPerFactor arrays.
 *
 * @param stExpression pointer to a list (length numFactors) of matrices
 * @param sparseRowsForFactor pointers into sparse structure for given factor (length numFactors + 1)
 * @param stMatrices length numFactors array of (double*) pointers to be filled with
 *        pointers to the contents of the matrices in ST.  Not used if NULL.
 * @param sparseMatrixStructure struct containing the dimensions of each factor, numGroupsPerFactor, 
 *        and maxFactorDimension
 */
void
getSparseContentAndStructure(const SEXP stExpression, const int *sparseRowForFactor,
                             double **stMatrices, SparseMatrixStructure *sparseMatrixStructure)
{
  int maxFactorDimension = 0;
  int numFactors = LENGTH(stExpression);
  
  for (int factor = 0; factor < numFactors; ++factor) {
    SEXP st_i = VECTOR_ELT(stExpression, factor);
    int factorDimension = INTEGER(getAttrib(st_i, R_DimSymbol))[0];
    
    if (factorDimension > maxFactorDimension) maxFactorDimension = factorDimension;
    
    if (stMatrices != NULL) stMatrices[factor] = REAL(st_i);
    
    sparseMatrixStructure->factorDimensions[factor] = factorDimension;
    
    int numRowsForFactor = sparseRowForFactor[factor + 1] - sparseRowForFactor[factor];
    sparseMatrixStructure->numGroupsPerFactor[factor] = numRowsForFactor / factorDimension;
  }
  
  sparseMatrixStructure->maxFactorDimension = maxFactorDimension;
}

/**
 * Return the index of the factor associated with parameter index ind
 *
 * @param row an index in [0, sparseRowsForFactor[numFactors] - 1]
 * @param numFactors total number of grouping factors
 * @param sparseRowsForFactor pointers into sparse structure for given factor (length numFactors + 1)
 *
 * @return index of assicated factor
 */
int getFactorForSparseRow(int row, int numFactors, const int *sparseRowsForFactor)
{
  
  for (int factor = 0; factor < numFactors; ++factor) {
    if (row < sparseRowsForFactor[factor + 1]) return factor;
  }
  error("invalid row index %d (max is %d)", row, sparseRowsForFactor[numFactors]);
  return -1;
}


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
applyPermutation(double *destination, const double *source, const int *permutation, int numValues)
{
  if (permutation == NULL) {
    Memcpy(destination, source, numValues);
  } else {
    for (int i = 0; i < numValues; ++i) {
      destination[i] = source[permutation[i]];
    }
  }
  return destination;
}

/**
 * Create PAX in dest.
 *
 * @param dest values to be calculated
 * @param perm NULL or a 0-based permutation vector defining P
 * @param A sparse design matrix
 * @param X dense design matrix
 * @param nc number of columns in X
 *
 */
void
multiplyWithPermutation(double *destination, const int *permutation,
                        const CHM_SP sparseModelMatrix,
                        const double *denseDesignMatrix, int numDenseColumns)
{
  int    *nonZeroRowIndices = (int *)    (sparseModelMatrix->i);
  int    *indicesForColumn  = (int *)    (sparseModelMatrix->p);
  double *sparseValues      = (double *) (sparseModelMatrix->x);
  
  // definitely worth highlighting right now that
  // the sparse model matrix is transposed, so that if it is n x q,
  // the transpose is q x n and its numCols is the same as the 
  // num rows in the dense matrix
  int numSparseRows = sparseModelMatrix->nrow;
  int numDenseRows  = sparseModelMatrix->ncol;
  
  // compute columns of the product as that aligns with how
  // the sparse matrix is stored
  double *columnOfProduct = Alloca(numSparseRows, double);
  R_CheckStack();
  
  for (int denseCol = 0; denseCol < numDenseColumns; ++denseCol) { // as many as num fixed effects
    AZERO(columnOfProduct, numSparseRows);
    
    for (int denseRow = 0; denseRow < numDenseRows; ++denseRow) { // as many as num observations
	    for (int sparseIndex = indicesForColumn[denseRow]; sparseIndex < indicesForColumn[denseRow + 1]; ++sparseIndex) {
      
        columnOfProduct[nonZeroRowIndices[sparseIndex]] +=
          denseDesignMatrix[denseRow + denseCol * numDenseRows] * sparseValues[sparseIndex];
      }
    }
    applyPermutation(destination + denseCol * numSparseRows,
                     columnOfProduct,
                     permutation,
                     numSparseRows);
  }
}

/**
 * Update the sparse model matrix A from the Zt and ST slots, where
 * A = S'T'Z'.
 *
 * @param regression an mer object
 */
void rotateSparseDesignMatrix(SEXP regression)
{
  CHM_SP sparseDesignMatrix        = Zt_SLOT(regression);
  CHM_SP rotatedSparseDesignMatrix = A_SLOT(regression);
  
  int numFactors = DIMS_SLOT(regression)[nt_POS];
  
  int *factorDimensions   = Alloca(numFactors, int);
  int *numGroupsPerFactor = Alloca(numFactors, int);
  double **stMatrices     = Alloca(numFactors, double *);
  R_CheckStack();
  
  int *sparseRowForFactor = Gp_SLOT(regression);
  
  int *indicesIntoZForColumn = (int*) (sparseDesignMatrix->p);
  int *indicesIntoAForColumn = (int*) (rotatedSparseDesignMatrix->p);
  
	int *nonZeroRowIndicesOfZ  = (int*) (sparseDesignMatrix->i);
  int *nonZeroRowIndicesOfA  = (int*) (rotatedSparseDesignMatrix->i);
  
  int numNonZeroesInZ = indicesIntoZForColumn[sparseDesignMatrix->ncol];
  int numNonZeroesInA = indicesIntoAForColumn[rotatedSparseDesignMatrix->ncol];
  
  double *valuesOfZ = (double *) (sparseDesignMatrix->x);
  double *valuesOfA = (double *) (rotatedSparseDesignMatrix->x);
  
  
  SparseMatrixStructure sparseMatrixStructure = { factorDimensions, numGroupsPerFactor, 0 };
  
  getSparseContentAndStructure(GET_SLOT(regression, lme4_STSym),
                               sparseRowForFactor, stMatrices,
                               &sparseMatrixStructure);
  
  // Copy Z' to A unless A has new nonzeros
  if (numNonZeroesInA == numNonZeroesInZ) {
    Memcpy(valuesOfA, valuesOfZ, numNonZeroesInZ);
  } else { // Only for nonlinear models with correlated random effects
    AZERO(valuesOfA, numNonZeroesInA); 	// Initialize potential nonzeros to 0
    
    for (int col = 0; col < rotatedSparseDesignMatrix->ncol; ++col) { /* Iterate over columns */
	    int aIndex = indicesIntoAForColumn[col];
      
	    for (int zIndex = indicesIntoZForColumn[col]; zIndex < indicesIntoZForColumn[col + 1]; ++zIndex) { // nonzeros in Z'
        while (nonZeroRowIndicesOfA[aIndex] < nonZeroRowIndicesOfZ[zIndex]) ++aIndex;          // matching pos in A
        if (nonZeroRowIndicesOfA[aIndex] != nonZeroRowIndicesOfZ[zIndex])
          error("nonconforming Zt and A structures, j = %d", col);
        valuesOfA[aIndex] = valuesOfZ[zIndex];
	    }
    }
  }
  /* When T != I multiply A on the left by T' */
  if (sparseMatrixStructure.maxFactorDimension > 1) {
    blockwiseLeftMultiplyMatrixByUpperUnitTrianglularMatrices(rotatedSparseDesignMatrix,
                                                              sparseRowForFactor,
                                                              stMatrices,
                                                              numFactors,
                                                              &sparseMatrixStructure);
  }
  blockwiseLeftMultiplyMatrixByDiagonalMatrices(rotatedSparseDesignMatrix,
                                                sparseRowForFactor, 
                                                stMatrices,
                                                numFactors,
                                                &sparseMatrixStructure);
}

/**
 * Multiply A on the left by T'
 *
 * @param A sparse model matrix, initialized to Zt
 * @param Gp group pointers
 * @param nc number of columns per term
 * @param nlev number of factors per term
 * @param st ST arrays for each term
 * @param nt number of terms
 * 
 */
static void
blockwiseLeftMultiplyMatrixByUpperUnitTrianglularMatrices(CHM_SP sparseDesignMatrix,
                                                          int *sparseRowForFactor, 
                                                          double **stMatrices, int numFactors,
                                                          const SparseMatrixStructure *sparseMatrixStructure)
{
  int *factorDimensions   = sparseMatrixStructure->factorDimensions;
  int *numGroupsPerFactor = sparseMatrixStructure->numGroupsPerFactor;
  
  int *nonZeroRowIndices = (int *)(sparseDesignMatrix->i);
  int *indicesForColumn  = (int *)(sparseDesignMatrix->p);
  double *values = (double *)(sparseDesignMatrix->x);
  
  double d_one[] = { 1.0, 0 };
  
  for (int col = 0; col < sparseDesignMatrix->ncol; ++col) {
    // multiply column j by T'
    for (int sparseIndex = indicesForColumn[col]; sparseIndex < indicesForColumn[col + 1];) {
	    int factorIndex = getFactorForSparseRow(nonZeroRowIndices[sparseIndex], numFactors, sparseRowForFactor);
	    
	    if (factorDimensions[factorIndex] <= 1) ++sparseIndex; // "T" is ident if dim is 1
	    else {
        // nudge this over until it's past this particular factor
        int maxSparseIndex = sparseIndex;
        while ((nonZeroRowIndices[maxSparseIndex] - sparseRowForFactor[factorIndex]) < numGroupsPerFactor[factorIndex]) {
          ++maxSparseIndex;
        }
        
        // number of rows in `B' in dtrmm call
        int numRows = maxSparseIndex - sparseIndex;	// numRows == 1 except in models with carry-over
        
        
        // perform B = B * T
        // the "U" prevents multiplication by S, as it considers ST as unit triangular
        F77_CALL(dtrmm)("R", "L", "N", "U",
                        &numRows,
                        factorDimensions + factorIndex,
                        d_one, stMatrices[factorIndex],
                        factorDimensions + factorIndex,
                        values + sparseIndex,
                        &numRows);
        sparseIndex += (numRows * factorDimensions[factorIndex]);
	    }
    }
  }
}


static void
blockwiseLeftMultiplyMatrixByDiagonalMatrices(CHM_SP sparseDesignMatrix,
                                              int *sparseRowForFactor,
                                              double **stMatrices, int numFactors,
                                              const SparseMatrixStructure *sparseMatrixStructure)
{
  int *factorDimensions   = sparseMatrixStructure->factorDimensions;
  int *numGroupsPerFactor = sparseMatrixStructure->numGroupsPerFactor;
  
  int *nonZeroRowIndices = (int *) (sparseDesignMatrix->i);
  int *indicesForColumn  = (int *) (sparseDesignMatrix->p);
  double *values = (double *) (sparseDesignMatrix->x);
  
  int numNonZeroes = indicesForColumn[sparseDesignMatrix->ncol];
    
  // Multiply A on the left by S
  for (int sparseIndex = 0; sparseIndex < numNonZeroes; ++sparseIndex) {
    // for each non-zero in A, find the factor to which it corresponds
    int factorIndex = getFactorForSparseRow(nonZeroRowIndices[sparseIndex], numFactors, sparseRowForFactor);
    
    // for that factor, find which diagonal element of S we should use in multiplication
    int stIndex = ((nonZeroRowIndices[sparseIndex] - sparseRowForFactor[factorIndex]) / numGroupsPerFactor[factorIndex]) *
      (factorDimensions[factorIndex] + 1);
    
    values[sparseIndex] *= stMatrices[factorIndex][stIndex];
  }
}
