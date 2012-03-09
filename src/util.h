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

// #define PRINT_TRACE
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

// functions that manipulate SEXPs

// alloc slot allocates a vector of the given type and length, stores it in
// the object's slot, and returns it
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length);
SEXP SET_DIMS(SEXP obj, int numRows, int numCols);
SEXP getListElement(SEXP list, char const* name);

#endif // BLME_UTIL_H
