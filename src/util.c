#include "util.h"

#include <Rmath.h> // for fabs

int allApproximatelyEqual(const double *p1, const double *p2, int numParameters, double tolerance)
{
  for (const double *p1End = p1 + numParameters; p1 < p1End; ) {
    if (fabs(*p1++ - *p2++) > tolerance) return 0;
  }
  
  return 1;
}

int allApproximatelyAbsolutelyEqual(const double *p1, const double *p2, int numParameters, double tolerance)
{
  for (const double *p1End = p1 + numParameters; p1 < p1End; ) {
    if (fabs(fabs(*p1++) - fabs(*p2++)) > tolerance) return 0;
  }
  
  return 1;
}

int allEqual(const double *p1, const double *p2, int numParameters)
{
  for (const double *p1End = p1 + numParameters; p1 < p1End; ) {
    if (*p1++ != *p2++) return 0;
  }
  
  return 1;
}

void printMatrix(const double *matrix, int numRows, int numCols)
{
  for (int row = 0; row < numRows; ++row) {
    for (int col = 0; col < numCols; ++col) {
      Rprintf("%f%s", matrix[row + col * numRows], (col < numCols ? " " : ""));
    }
    Rprintf("\n");
  }
}
