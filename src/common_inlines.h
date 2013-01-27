// DO NOT ADD ifdef inclusion rules here.
// Also, do not include in any header file.

// This file should be included directly into C files,
// as it should contain only *static* functions and variables.
/**
 * Return the sum of squares of the first n elements of x
 *
 * @param n
 * @param x
 *
 * @return sum of squares
 */
static R_INLINE
double getSumOfSquares(const double *x, int n) {
  double result = 0.0;
  for (int i = 0; i < n; ++i) result += x[i] * x[i];
  return(result);
}
