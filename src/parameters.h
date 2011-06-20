/* supports the conversion between different parameter types
 */

#ifndef BLME_PARAMETERS_H
#define BLME_PARAMETERS_H

typedef enum {
  PARAMETERIZATION_COVARIANCE,     // covariance is just the matrix
  PARAMETERIZATION_SD_CORRELATION, // [ sds, lower triangle of correlation ]
  PARAMETERIZATION_SPECTRAL,       // [ eigen values, angles for givens rotations ]
  PARAMETERIZATION_ST,             // [ S, lower tri of T ]
  PARAMETERIZATION_CHOLESKY,       // cholesky is simply the result of a call to LAPACK for a decomp; is matrix
  PARAMETERIZATION_PRIOR,          // depends on prior
} parameterizationType_t;

/* here is a detailed definition of all of the different parameterizations
 * 
 * COVARIANCE:
 *   length  : dim^2
 *   contents: actual covariance matrix
 *
 * SD_CORRELATION:
 *   length  : dim + dim * (dim - 1) / 2
 *   contents: first dim components are the standard deviations
 *             second (dim choose 2) components are the lower triangle of the correlation matrix
 *
 * SPECTRAL:
 *   length  : dim + dim * (dim - 1) / 2
 *   contents: first dim components are the eigenvalues
 *             second (dim choose 2) components is the lower triangle of the skew-symmetric matrix
 *             of a Cayley transformation
 *
 * ST:
 *   length  : dim + dim * (dim - 1) / 2
 *   contents: first dim components are diagonal of S
 *             second (dim choose 2) components are the lower triangle of unit triangular matrix T
 *             covariance matrix is TSST'
 *
 * CHOLESKY:
 *   length  : dim^2
 *   contents: column major matrix containing the lower triangle of the cholesky
 *             decomposition, result of a call to F77_CALL(dpotrf)
 *
 * PRIOR (spectral):
 *   contents: same as spectral
 *
 * PRIOR (correlation):
 *   length  : dim + dim + dim * (dim - 1) / 2
 *   contents: first dim components are the elements of the scale (diagonal) matrix
 *             remaining components are an ST decomp of the inner matrix
 */
 
/*
 * Some of these conversions are implemented directly, some go through other representations.
 * Those that are direct are marked.
 * Typically, those that are not go through the covariance matrix parameterization. For
 * some of the indirect implementations, it may be possible to increase the speed by hand-
 * rolling the implementation.
 */
  
/*
 * Finally, those that take (double *) as a source instead of (const double *) use the source
 * as scratch.
 */

void convertSDCorrelationToCovariance(const double *source, int dim, double *target);       // direct
void convertSpectralToCovariance(const double *source, int dim, double *target);            // direct
void convertSTToCovariance(const double *source, int dim, double *target);                  // direct
void convertSTToCovarianceWithScale(const double *source, int dim, double *target,          // direct
                                    double matrixScale);
void convertCholeskyToCovariance(const double *source, int dim, double *target);            // direct
void convertPriorCorrelationToCovariance(const double *source, int dim, double *target);    // direct

void convertSTToCovarianceInverse(const double *source, int dim, double *target);           // direct
void convertSTToCovarianceInverseWithScale(const double *source, int dim, double *target,   // direct
                                           double matrixScale);


void convertCovarianceToSDCorrelation(const double *source, int dim, double *target);       // direct
void convertSpectralToSDCorrelation(const double *source, int dim, double *target);
void convertSTToSDCorrelation(const double *source, int dim, double *target);
void convertCholeskyToSDCorrelation(const double *source, int dim, double *target);
void convertPriorCorrelationToSDCorrelation(const double *source, int dim, double *target); // direct

void convertCovarianceToSpectral(double *source, int dim, double *target);                  // direct, destroys source
void convertSDCorrelationToSpectral(const double *source, int dim, double *target);
void convertSTToSpectral(const double *source, int dim, double *target);
void convertCholeskyToSpectral(const double *source, int dim, double *target);
void convertPriorCorrelationToSpectral(const double *source, int dim, double *target);

void convertCovarianceToST(double *source, int dim, double *target);                        // destroys source
void convertSDCorrelationToST(const double *source, int dim, double *target);
void convertSpectralToST(const double *source, int dim, double *target);
void convertCholeskyToST(const double *source, int dim, double *target);                    // direct
void convertPriorCorrelationToST(const double *source, int dim, double *target);            // direct

void convertCovarianceToCholesky(const double *source, int dim, double *target);            // direct
void convertSDCorrelationToCholesky(const double *source, int dim, double *target);
void convertSpectralToCholesky(const double *source, int dim, double *target);
void convertSTToCholesky(const double *source, int dim, double *target);                    // direct
void convertPriorCorrelationToCholesky(const double *source, int dim, double *target);

void convertCovarianceToPriorCorrelation(double *source, int dim, double *target);          // destroys source
void convertSDCorrelationToPriorCorrelation(const double *source, int dim, double *target);
void convertSpectralToPriorCorrelation(const double *source, int dim, double *target);
void convertSTToPriorCorrelation(const double *source, int dim, double *target);            // direct
void convertCholeskyToPriorCorrelation(const double *source, int dim, double *target);

#endif /* BLME_PARAMETERS_H */
