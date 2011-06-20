#ifndef BLME_UTIL_H
#define BLME_UTIL_H

// sadly, this needs to be here before any redefinitions of alloca,
// as it includes alloca itself via <stdlib.h>
#include <R.h>
#include <Rdefines.h>

/* When appropriate, alloca is cleaner than malloc/free.  The storage
 * is freed automatically on return from a function. When using gcc the
 * builtin version is much faster. */
#ifdef __GNUC__
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
/* this is necessary (and sufficient) for Solaris 10: */
# ifdef __sun
#  include <alloca.h>
# endif
#endif

#ifdef PRINT_TRACE
#  define DEBUG_PRINT_ARRAY(header, array, length) { \
     unsigned long long *_X_ = (unsigned long long *) (array); \
     int _I_, _SZ_ = (length); \
     Rprintf("%s: %llu", (header), _X_[0]); \
     for(_I_ = 1; _I_ < _SZ_; ++_I_) Rprintf(" %llu", _X_[_I_]); \
     Rprintf("\n"); \
   }
#else
#  define DEBUG_PRINT_ARRAY(header, array, length)
#endif

/** zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/** alloca n elements of type t */
#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

int allApproximatelyEqual(const double *p1, const double *p2, int numParameters, double tolerance);
int allApproximatelyAbsolutelyEqual(const double *p1, const double *p2, int numParameters, double tolerance);
int allEqual(const double *p1, const double *p2, int length);

void printMatrix(const double *matrix, int numRows, int numColumns);

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
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
  SEXP val = allocVector(type, length);
  
  SET_SLOT(obj, nm, val);
  return val;
}

static R_INLINE
SEXP SET_DIMS(SEXP obj, int numRows, int numCols)
{
  SEXP dimsExp = allocVector(INTSXP, 2);
  int *dims = INTEGER(dimsExp);
  dims[0] = numRows;
  dims[1] = numCols;
  
  setAttrib(obj, R_DimSymbol, dimsExp);
  
  return(obj);
}

#endif // BLME_UTIL_H
