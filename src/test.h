#ifndef BLME_TEST_H
#define BLME_TEST_H

#include <R.h>
#include <Rinternals.h>

SEXP bmer_matrixTest();
SEXP bmer_parametersTest();
SEXP bmer_lmmTest();

SEXP createTestRegression(); // defined in lmer_test.c

#endif // BLME_TEST_H
