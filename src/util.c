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

// ripped from Matrix package
/**
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.  The validity of this function
 * depends on SET_SLOT not duplicating val when NAMED(val) == 0.  If
 * this behavior changes then ALLOC_SLOT must use SET_SLOT followed by
 * GET_SLOT to ensure that the value returned is indeed the SEXP in
 * the slot.
 * NOTE:  GET_SLOT(x, what)        :== R_do_slot       (x, what)
 * ----   SET_SLOT(x, what, value) :== R_do_slot_assign(x, what, value)
 * and the R_do_slot* are in src/main/attrib.c
 *
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 *
 * @return SEXP of given type and length assigned as slot nm in obj
 */
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
  SEXP val = allocVector(type, length);
  
  SET_SLOT(obj, nm, val);
  return val;
}

SEXP SET_DIMS(SEXP obj, int numRows, int numCols)
{
  SEXP dimsExp = allocVector(INTSXP, 2);
  int *dims = INTEGER(dimsExp);
  dims[0] = numRows;
  dims[1] = numCols;
  
  setAttrib(obj, R_DimSymbol, dimsExp);
  
  return(obj);
}

SEXP getListElement(SEXP list, char const* name)
{
  SEXP result = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    
  for (unsigned int i = 0; i < length(list); ++i)
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      result = VECTOR_ELT(list, i);
      break;
    }
  return result;
}
