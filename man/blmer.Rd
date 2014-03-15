\name{blme}
\alias{blmer}
\alias{bglmer}
\concept{GLMM}% << so it's found
\concept{NLMM}% << so it's found
\title{Fit Bayesian Linear and Generalized Linear Mixed-Effects Models}
\description{
  Maximum a posteriori estimation for linear and generalized
  linear mixed-effects models in a Bayesian setting. Built off of
  \code{\link[lme4]{lmer}}.
}
\usage{
blmer(formula, data, REML = TRUE,
      control = lmerControl(), start = NULL, verbose = 0L,
      subset, weights, na.action, offset, contrasts = NULL,
      devFunOnly = FALSE, cov.prior = wishart,
      fixef.prior = NULL, resid.prior = NULL, \dots)
bglmer(formula, data, family = gaussian,
       control = glmerControl(), start = NULL, verbose = 0L,
       nAGQ = 1L, subset, weights, na.action, offset,
       contrasts = NULL, mustart, etastart,
       devFunOnly = FALSE, cov.prior = wishart,
       fixef.prior = NULL, \dots)
}
\arguments{
  \item{cov.prior}{a BLME \link[=bmerDist-class]{prior} or list of priors with allowable
    distributions: \code{wishart}, \code{invwishart}, \code{gamma},
    \code{invgamma}, or \code{NULL}. Imposes a prior over the covariance of the random
    effects/modeled coefficients. Default is \code{wishart}. The \code{NULL} argument
    imposes flat priors over all relevant parameters.}
  \item{fixef.prior}{a BLME prior of family \code{normal}, \code{t}, or \code{NULL}.
    Imposes a prior over the fixed effects/modeled coefficients.
    Default is \code{NULL}.}
  \item{resid.prior}{a BLME prior of family \code{gamma}, \code{invamma}, \code{point}
    or \code{NULL}. Imposes a prior over the noise/residual variance, also known as common scale
    parameter or the conditional variance given the random effects.
    Default is \code{NULL}.}
  \item{start}{like the \code{start} arguments for \code{\link[lme4]{lmer}} and
    \code{\link[lme4]{glmer}} a numeric vector or named list. Unlike the aforementioned,
    list members of \code{fixef} and \code{sigma} are applicable to linear mixed models
    provided that numeric optimization is required for these parameters.}
  \item{formula, data, REML, family, control, verbose, nAGQ,
    mustart, etastart, devFunOnly, \dots}{model specification arguments as in \code{\link[lme4]{lmer}} and \code{\link[lme4]{glmer}};
    see there for details.}
  \item{subset, weights, na.action, offset, contrasts}{further model
    specification arguments as in \code{\link[stats]{lm}}; see there for
    details.}
}
\details{
  The bulk of the usage for \code{blmer} and \code{bglmer} closely
  follows the functions \code{\link[lme4]{lmer}} and
  \code{\link[lme4]{glmer}}. Those help pages provide a good overview of
  fitting linear and generalized linear mixed models. The primary
  distinction is that \code{blmer} and \code{bglmer} allow the user to
  do Bayesian inference or penalized maximum likelihood, with priors imposed on the different
  model components. For the specifics of any distribution listed below,
  see the \link[=bmerDist-class]{distributions} page.
  
  \strong{Covariance Prior}
  
  The \code{cov.prior} argument applies a prior over the
  covariance matrix of the random effects/modeled coefficients.
  As there is one covariance matrix for every named grouping factor -
  that is every element that appears to the right of a vertical bar
  ("|") in the model formula - it is possible to apply as many
  different priors as there are said factors.
  
  The general formats of an argument to \code{blmer} or \code{bglmer}
  for such a prior are of the form:

\itemize{
  \item \code{cov.prior = factor.name ~ covariance.distribution(option1 = value1, \dots)}
  \item \code{cov.prior = list(fc.nm ~ dist1, fc.nm ~ dist2, ..., default.distribution)}
}
  
  If the \dQuote{\code{factor.name ~}} construct is ommitted, the prior
  is interpretted as a default and applied to all factors that
  lack specific priors of their own. Options are not required,
  but permit fine-tuning of the model.

  Supported distributions are \code{gamma}, \code{invgamma}, \code{wishart},
  \code{invwishart}, \code{NULL}, and \code{custom}.

  The \code{common.scale} option, a logical, determines whether or
  not the prior applies to in the absolute-real world
  sense (value = \code{FALSE}), or if the prior is applied to the random effect
  covariance divided by the estimated residual variance (\code{TRUE}). As a practical matter,
  when false computation can be slower as the profiled common scale may
  no longer have a closed-form solution. As such, the default for all
  cases is \code{TRUE}.

  Other options are specified along with the specific distributions and
  defaults are explained in the blme \link[=bmerDist-class]{distributions} page.  

  \strong{Fixed Effects Prior}
  
  Priors on the fixed effects, or unmodeled coefficients, are specified
  in a fashion similar to that of covariance priors. The general format is

\itemize{\item \code{fixef.prior = multivariate.distribution(options1 = value1, \dots)}}
  
  At present, the implemented multivariate distributions are \code{normal}, \code{t},
  and \code{NULL}. \code{t} priors cannot be used when \code{REML} is \code{TRUE},
  as that integral does not have a closed form solution.

  \strong{Residual Variance Prior}

  The general format for a residual variance prior is the same as for a fixed
  effect prior. The supported distributions are \code{point}, \code{gamma},
  \code{invgamma}.

}
\value{
  An object of class \code{"\linkS4class{bmerMod}"}, for which many methods
  are available. See there for details.
}
\seealso{
  \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
  \code{\linkS4class{merMod}} class, and \code{\link[stats]{lm}}.
}
\examples{
##  covariance prior
(fm1 <- blmer(Reaction ~ Days + (0 + Days|Subject), sleepstudy,
              cov.prior = gamma))
(fm2 <- blmer(Reaction ~ Days + (0 + Days|Subject), sleepstudy,
              cov.prior = gamma(shape = 2, rate = 0.5, posterior.scale = 'sd')))
(fm3 <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
              cov.prior = wishart))
(fm4 <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
              cov.prior = invwishart(df = 5, scale = diag(0.5, 2))))

##  custom prior
penaltyFn <- function(sigma)
  dcauchy(sigma, 0, 10, log = TRUE)
(fm5 <- blmer(Reaction ~ Days + (0 + Days|Subject), sleepstudy,
              cov.prior = custom(penaltyFn, chol = TRUE, scale = "log")))


##  fixed effect prior
(fm6 <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
              cov.prior = NULL,
              fixef.prior = normal))
(fm7 <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
              cov.prior = NULL,
              fixef.prior = normal(cov = diag(0.5, 2), common.scale = FALSE)))

##  residual variance prior
##  eight schools example
y <- c(28, 8, -3, 7, -1, 1, 18, 12);
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18);
g <- 1:8;
  
(schools <- blmer(y ~ 1 + (1 | g), resid.prior = point,
                  cov.prior = NULL, REML = FALSE,
                  weights = 1 / sigma^2));
}
\keyword{models}
\keyword{methods}
