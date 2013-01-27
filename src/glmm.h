#ifndef BLME_GLMM_H
#define BLME_GLMM_H

#include <R.h>
#include <Rinternals.h>

#include "cache.h"

MERCache* createGLMMCache(SEXP regression);
void deleteGLMMCache(MERCache *cache);

#endif // BLME_GLMM_H
