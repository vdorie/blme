\name{blme}
\alias{blmer}
\alias{bglmer}
\concept{GLMM}% << so it's found
\concept{NLMM}% << so it's found
\title{Fit Bayesian Linear and generalized linear Mixed-Effects models}
\description{
  Maximum a posteriori estimation for linear and generalized
  linear mixed-effects models in a Bayesian setting. Extension
  of, and similar to, \code{\link[lme4]{lmer}}.
}
\usage{
blmer(formula, data, family = NULL, REML = TRUE,
      control = list(), start = NULL, verbose = FALSE,
      doFit = TRUE, subset, weights, na.action, offset,
      contrasts = NULL, model = TRUE, x = TRUE,
      cov.prior = "wishart", fixef.prior = "normal",
      var.prior = "inverse.gamma",
      \dots)
bglmer(formula, data, family = gaussian, start = NULL,
       verbose = FALSE, nAGQ = 1, doFit = TRUE, subset, weights,
       na.action, offset, contrasts = NULL, model = TRUE,
       control = list(), cov.prior = "wishart",
       fixef.prior = "normal", \dots)
}
\arguments{
  \item{cov.prior}{character string or list of strings. Imposes a prior over the
    covariance of the modeled coefficients. Currently supported
    types are \code{"none"}, \code{"spectral"}, \code{"correlation"},
    and the direct application of a named distribution, i.e.
    \code{"gamma"}, \code{"inverse.gamma"}, \code{"wishart"}, or
    \code{"inverse.wishart"}. Default is \code{"wishart"} with
    parameters as specified in the details section. The \code{NULL} argument
    imposes flat priors over all parameters.}
  \item{fixef.prior}{string. Imposes a prior over the modeled
    coefficients, also known as the fixed effects. Only currently
    supported options are \code{"none"} and \code{"normal"}.
    Default is \code{"normal"} with parameters as specified in the
    details section. The \code{NULL} argument imposes flat priors over all
    parameters.}
  \item{var.prior}{string. Imposes a prior over the so-called common scale
    parameter, also known as the conditional variance given the modeled
    coefficients/random effects. Options are \code{"point"}, \code{"gamma"}, and
    \code{"inverse.gamma"}. Default is \code{"inverse.gamma"} with parameters
    specified below.  The \code{NULL} argument imposes a flat prior.}
  \item{verbose}{logical scalar. If \code{TRUE} verbose output is
    generated during the optimization of the parameter estimates. In
    addition, debugging information about priors is printed out before
    optimization.}
  \item{formula, data, family, REML, nAGQ, control, start, doFit,
    model, x, ...}{model specification arguments as in
    \code{\link[lme4]{lmer}}; see there for details.}
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
  do Bayesian inference, with priors imposed on the different
  model components.
  
  \strong{Covariance Prior}
  
  The \code{cov.prior} argument applies a prior over the
  covariance matrix of the modeled coefficients, or "random effects".
  As there is one covariance matrix for every named grouping factor -
  that is every element that appears to the right of a vertical bar
  ("|") in the model formula - it is possible to apply as many
  different priors as there are said factors.
  
  The general format of an argument to \code{blmer} or \code{bglmer}
  for such a prior is of the form:
  
  \preformatted{
    cov.prior = "factor.name ~ prior.type(type.options)"
  }
  
  If the "\code{factor.name ~}" construct is ommitted, the prior
  is interpretted as a default and applied to all factors that
  lack specific priors of their own. Options are not required,
  but permit fine-tuning of the model.
  
  \code{prior.type} can be \code{none}, \code{correlation},
  \code{spectral}, or a direct prior of a named distribution.
  Respectively, these signify:
  \itemize{
    \item \code{none} or \code{NULL} - no/flat prior, useful to override a default
    \item \code{correlation} - decomposes each covariance matrix
      into a diagonal component and a center component, roughly
      corresponding to the standard deviations and correlations
      of the original matrix. If \eqn{D} is diagonal, \eqn{R} is
      a correlation matrix, then the decomposition is
      \eqn{\Sigma = DRD}. As estimation in space is highly constrained,
      \eqn{R} is typically taken to be any positive definite matrix.
    \item \code{spectral} - applies priors to a spectral/eigen decomposition
      of each covariance matrix, typically with a flat prior over
      the orthonormal component and idendependent and identical priors
      over the eigenvalues.
    \item direct - any of \code{gamma}, \code{inverse.gamma},
      \code{wishart}, or \code{inverse.wishart} family of distributions,
      to be applied directly to the covariance matrix, or possibly the
      standard deviation in the univiariate case.
  }

  Type-specific options include the further detailing of the distributions
  to be used in a decomposition, as well as a few more general settings.
  An option for all of cases is
  \code{common.scale}, which accepts a boolean argument or the strings
  \code{'true'} or \code{'false'}. When true, the prior is assumed to
  be over the covariance of the random effects with the marginal scale
  component factored out. Conversely, when false, the prior is applied
  over the covariance in an absolute sense, as would be appropriate
  for substantive prior knowledge. Practically, when false computation is
  more difficult as profiling the common scale may no longer have
  a closed-form solution. As such, the default for all cases is \eqn{TRUE}.

  The \code{data.scale} setting can be either \code{'absolute'} or \code{'free'}.
  When set to absolute, scale (and inverse scale/rate) parameters are interpretted
  as they are given. Using \code{data.scale = 'free'} adjusts the scales
  of each prior by \code{sd(y) / sd(x)} (or the same with variances, if applicable),
  so that they correspond to an absolute prior specified on standardized data.
  Exceptions to this rule are that \code{sd(y)} is replaced by \eqn{1} when the
  model is not over continuous outcomes, and the same change follows for input
  variables that have no variance, like intercepts. The quantity
  \code{sd(y) / sd(x)}, with suitable substitions, will be referred to as the
  \code{sd.ratio} hereafter.

  Further options for \code{correlation} and \code{spectral} types are the
  distributional assumptions that are necessary for them to be fully specified.
  By default the \code{correlation} type applies independent \code{gamma} priors to each
  component of \eqn{D}, and a \code{wishart} prior to \eqn{R}. The default
  for \code{spectral} is to apply independent \code{gamma} priors over the
  eigenvalues.
  
  To change these settings, these prior types take options of the form:
  
  \preformatted{
    coordinate.name ~ univariate.distribution(distribution.options)
  }
  
  In the case of the correlation type, also possible is:
  
  \preformatted{
    multivariate.distribution(distribution.options)
  }
  
  Coordinate names (or numbers) refer to the components of the covariance
  matrix. For example, a model containing \code{(1 + x.1 | g.1)} has
  coordinates \code{(Intercept)} and \code{x.1} at grouping factor \code{g.1}.
  Eigenvalues do not directly correspond to coordinates as named, so
  they must be specified by the number of their position along the diagonal.
  
  Options for different distributions are specified in a named-list style,
  and are comprised of:
  
  \itemize{
    \item \code{gamma} - \code{shape}, \code{rate}, and \code{posterior.scale}.
    Shape and rate are nonnegative real numbers, while the posterior scale can be
    \code{'sd'} or \code{'var'}, specifying on which form of the parameter the
    prior is to be applied.
    \item \code{inverse.gamma} - \code{shape}, \code{scale}, and
    \code{posterior.scale}. Arguments are similar to the \code{gamma} case.
    \item \code{wishart} - \code{df}, \code{scale}. The scale argument must
    be a positive semidefinite matrix of the appropriate dimension, a nonnegative
    scalar, or vector of nonnegative scalars which is turned into a diagonal
    matrix.
    \item \code{inverse.wishart} - \code{df}, \code{inverse.scale}. Similar
    to the \code{wishart} case.
  }
  
  For the \code{correlation} type, family defaults are:
  \itemize{
    \item \code{gamma} - \code{shape = 2}, \code{posterior.scale = 'sd'}.
      A partial default of \code{shape = 1.5} when
      \code{posterior.scale = 'var'} is set by itself.
      \code{rate} defaults a value that puts the mode at \eqn{10^1} or
      \eqn{10^2}, for \code{posterior.scale = 'sd'} and \code{'var'}
      respectively. If the mode does not exist, the mean is used instead.
      Unlike the standard behavior, for \code{data.scale = 'free'},
      the square root of \code{sd.ratio} or \code{sd.ratio} itself further
      divides the rate. The square roots of
      a standard deviation is used as for an  unconstrained positive definite
      matrix for \eqn{R}, the full covariance matrix, \eqn{DRD}, will have a
      diagonal consisting of the squares of two numbers, so that "scales"
      will appear four times. 
    \item \code{inverse.gamma} - \code{shape = 0.5},
      \code{posterior.scale = 'sd'}, and the \code{scale} is chosen to set
      the mode equal \eqn{10^1} or \eqn{10^2}, depending on the
      value of \code{posterior.scale} as above.
    \item \code{wishart} - for a grouping factor of dimension \eqn{K}, 
      \code{df = K + 1}. When the mode of the distribution exists, the
      \code{scale} is chosen to set the mode equal \eqn{10^2} times the
      identity matrix. When, the mode does not exist, the mean is
      set to that value.
    \item \code{inverse.wishart} - for a grouping factor of dimension \eqn{K},
      \code{df = K - 0.5}. \code{inverse.scale} is chosen so that the mode is
      equal to \eqn{10^2} times the identity matrix.
  }
  
  For the \code{spectral} type:
  \itemize{
    \item \code{gamma} - \code{shape = 1.5}, \code{posterior.scale = 'var'},
      with a partial default of \code{shape = 2} when \code{posterior.scale = 'sd'}
      is specified by itself. \code{rate} defaults to setting the mode to \eqn{10^4}
      or \eqn{10^2}, depending on whether \code{posterior.scale = 'var'} or \code{'sd'},
      respectively. If the mode does not exist, the mean is set to that value instead.
    \item \code{inverse.gamma} - \code{shape = 0.5},
      \code{posterior.scale = 'var'}, and the \code{scale} is chosen to set
      the mode equal to \eqn{10^4} or \eqn{10^2},
      depending on the value of \code{posterior.scale}.
  }
  
  Finally, for directly applied distributions:
  \itemize{
    \item \code{gamma} - \code{shape = 2}, \code{posterior.scale = 'sd'},
      with a partial default of \code{shape = 1.5} when \code{posterior.scale = 'var'}
      is set by itself. \code{rate} is chosen so that the mode, when it exists,
      is equal to \eqn{10^2} (\eqn{10^4} for \code{posterior.scale = 'var'}). When
      the mode does not exist, the mean is fixed to that value instead.
    \item \code{inverse.gamma} - \code{shape = 0.5}, \code{posterior.scale = 'sd'},
      and the \code{scale} is chosen so that the mode is equal to \eqn{10^2} or
      \eqn{10^4}, for \code{posterior.scale = 'sd'} and \code{'var'}, respectively.
    \item \code{wishart} - for a grouping factor of dimension \eqn{K}, \code{df = 1 + 3}.
      \code{scale} is chosen so that the mode of the distribution is equal to
      \eqn{10^4} times the identity matrix, whenever it is that the mode exists.
      When it does not, the mean is set to that value instead.
    \item \code{inverse.wishart} - for a grouping factor of dimension \eqn{K},
      \code{df = K - 0.5}. \code{inverse.scale} is chosen to set the mode of the 
      distribution to \eqn{10^4} times the identity matrix.
  }
  
  \strong{Fixed Effects Prior}
  
  Priors on the 'fixed effects', or unmodeled coefficients, are specified
  in a fashion similar to that of covariance priors. The general format is
  
  \preformatted{
    fixef.prior = "multivariate.distribution(distribution.options)"
  }
  
  At present, the only implemented multivariate distributions are \code{normal}
  and \code{none}. The \code{normal} distribution takes as options:
  \itemize{
    \item \code{sd} - a single scalar, a vector of length 2, or a vector of
      length equal to the number of unmodeled coefficients. Provides the
      standard deviations for an assumed to be independent set of Gaussian priors.
      All elements must be nonnegative. When a single value is applied, it is
      repeated for each unmodeled coefficient. When two values are given, the
      first is used only for the intercept term, while the other is repeated. 
    \item \code{cov} - same as above, but also permits specifying full matrices
      of the appropriate dimension. Note that \code{sd} is on the scale of
      standard deviations, while \code{cov} is that of a variance/covariance.
      Must be positive semidefinite, or, for vector form, consistingly entirely of
      nonnegative numbers.
    \item \code{common.scale} - a boolean, as above.
      Optimization in \code{\link[lme4]{lmer}} is done with all covariances
      multiplied by a "common scale factor". When the \eqn{TRUE}, it is assumed that
      the specified covariance should also be multiplied by this factor, simplifying
      computation. \eqn{FALSE} permtis the specification of the covariance in an
      absolute, real-world sense.
    \item \code{data.scale} - either of \code{'absolute'} or \code{'free'}.
      The absolute setting interprets standard deviations and covariances
      as given, while a prior that is scale-free has a covariance matrix that
      is scaled by \code{sd.ratio} squared.
  }
  
  The default is for \code{sd} to be \code{c(10, 2.5)}, \code{common.scale} to be
  \eqn{TRUE}, and \code{data.scale} to be \code{'free'}.

  \strong{Common Scale Prior}

  The general format for a common scale prior is the same as for a direct prior
  above. Supported families with options are:

  \itemize{
    \item \code{none}
    \item \code{point} - fixes the common scale to a particular value - specified by
      name as \code{value}. The \code{posterior.scale} option, being either \code{'sd'}
      \code{'var'}, specifies the scale of \code{value}. With no value specified,
      defaults to \eqn{1}.
    \item \code{gamma} - Options are \code{shape}, \code{rate}, and \code{posterior.scale},
      with the same interpretations as for the covariance priors above.
    \item \code{inverse.gamma} - Options are \code{shape}, \code{scale}, and
      \code{posterior.scale}.
  }

The default is an improper inverse gamma distribution with both parameters as 0, i.e. proportional to \eqn{\sigma^-1}.

}
\value{
  An object of class \code{"\linkS4class{bmer}"}, for which many methods
  are available.  See there for details.
}
\seealso{
  \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
  \code{\linkS4class{mer}} class, and \code{\link[stats]{lm}}.
}
\examples{
#  covariance prior
(fm1 <- blmer(Reaction ~ Days + (Days|Subject), sleepstudy,
              cov.prior = "gamma"))
(fm2 <- blmer(Reaction ~ Days + (Days|Subject), sleepstudy,
              cov.prior = "gamma(shape = 2, rate = 0.5, posterior.scale = 'sd')"))
(fm3 <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
              cov.prior = "wishart"))
(fm4 <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
              cov.prior = "inverse.wishart(df = 5, inverse.scale = diag(0.5, 2))"))

#  unmodeled coefficient prior
(fm5 <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
              cov.prior = NULL, var.prior = NULL,
              fixef.prior = "normal"))
(fm6 <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
              cov.prior = NULL, var.prior = NULL,
              fixef.prior = "normal(cov = diag(0.5, 2), common.scale = FALSE)"))

#  common scale prior
#    eight schools example
  y <- c(28, 8, -3, 7, -1, 1, 18, 12);
  sigma <- c(15, 10, 16, 11, 9, 11, 10, 18);
  y.z <- (y - mean(y)) / sigma;
  g <- 1:8
  
  (schools <- blmer(y.z ~ 1 + (1 | g), var.prior = "point",
                   cov.prior = NULL, fixef.prior = NULL));
}
\keyword{models}
\keyword{methods}
