#ifndef BLME_DIST_GAMMA_H
#define BLME_DIST_GAMMA_H

// log-normalizing constants return 0.0 for improper dists
// values at 0 are defined below. for negative inputs, return NaN

// gammas are propto x^(shape - 1) exp(-rate * x)

// impropers are allowed, with the following rules:
//   shape >= 0, rate >= 0
//   shape > 1 => p(0) := 0

// inverse gammas are propto x^-(shape + 1) exp(-scale / x)
//
// improper rules:
//   shape >= 0, scale >= 0
//   scale == 0 => p(0) := inf

// normalizing constant
double gammaDevianceConstantPart(double shape, double rate);
// part that depends on x
double gammaDevianceVaryingPart(double x, double shape, double rate);

// these are literal, and don't consider the limiting behavior of their counterparts
double gammaDevianceVaryingPolynomialPart(double x, double shape);
double gammaDevianceVaryingExponentialPart(double x, double rate);


double inverseGammaDevianceConstantPart(double shape, double scale);
double inverseGammaDevianceVaryingPart(double x, double shape, double scale);
double inverseGammaDevianceVaryingPolynomialPart(double x, double shape);
double inverseGammaDevianceVaryingExponentialPart(double x, double scale);


#endif /* BLME_DIST_GAMMA_H */
