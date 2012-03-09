#ifndef BLME_SIM_H
#define BLME_SIM_H

#include <R.h>
#include <Rdefines.h>

// Simulates from an approximate posterior.
SEXP bmer_sim(SEXP regression, SEXP numSims);

#endif
