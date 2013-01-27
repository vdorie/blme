#ifndef BLME_DIST_MULTIVARIATE_NORMAL_H
#define BLME_DIST_MULTIVARIATE_NORMAL_H

/**
 * Density of a multivariate normal with given mean and covariance
 *
 * @param x            normally distributed vector for which the density is desired
 * @param dim          length of x, the mean, and the sides of the covar
 * @param mean         mean of x. can be NULL
 * @param sdInverse    1 / single standard deviation to be used for all (iid)
 * @param sdsInverse   1 / vector of sds (indep)
 * @param covInverse   arbitrary covariance matrix
 * @param logDetCov    the log of the determinant of the covariance matrix
 * @param useLogScale  whether or not the result is returned as the logarithm of the density
 *
 * @return the density or log density of x
 */
double dmvn(const double *x, int dim, const double *mean, double sdInverse,
            double logDetCov, int useLogScale);
double dmvn2(const double *x, int dim, const double *mean, const double* sdsInverse,
             double logDetCov, int useLogScale);
double dmvn3(const double *x, int dim, const double *mean, const double* covInverse,
             double logDetCov, int useLogScale);

#endif /* BLME_DIST_MULTIVARIATE_NORMAL_H */
