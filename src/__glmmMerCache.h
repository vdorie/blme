#ifndef __GLMM_MER_CACHE__
#define __GLMM_MER_CACHE__

// supposedly opaque class; not to be mucked around with

struct _MERCache {
  // anything that needs to be shared between lmm and glmm should come first
  // these *NEED* to be __glmmMerCache.h and __merCache.h
  double priorDevianceConstantPart;
};

#endif
