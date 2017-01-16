\docType{class}
\name{bmerDist-class}
\title{Bayesian Linear Mixed-Effects Model Prior Representations and bmer*Dist Methods}
% Classes
\alias{bmerDist}
\alias{bmerDist-class}
% Methods
\alias{print,bmerDist-method}
\alias{show,bmerDist-method}

\description{Objects created in the initialization step of a \pkg{blme}
  model that represent the type of prior being applied.}

\section{Objects from the Class}{Objects can be created by calls of the
  form \code{new("bmerPrior", ...)}  or, more commonly, as side effects of the
  \code{\link{blmer}} and \code{\link{bglmer}} functions.

  When using the main \code{blme} functions, the prior-related arguments can be
  passed what essentially
  are function calls with the distinction that they are delayed in evaluation
  until information about the model is available. At that time, the functions
  are \emph{defined} in a special environment and then \emph{evaluated} in an
  environment that directly inherits form the one in which \code{blmer} or
  \code{bglmer} was called. This is reflected in some of the
  prototypes of various prior-creating functions which depend on parameters not
  available in the top-level environment.

  Finally, if the trailing parentheses are omitted from a \code{blmer}/\code{bglmer}
  prior argument, they are simply added as a form of \dQuote{syntactic sugar}.
}

\section{Prior Distributions}{
  This section lists the prototypes for the functions that are called to parse
  a prior during a model fit.

  \strong{Fixed Effect Priors}

  \itemize{
    \item \code{normal(sd = c(10, 2.5), cov, common.scale = TRUE)}
    
    Applies a Gaussian prior to the fixed effects. Normal priors are constrained
    to have a mean of 0 - non-zero priors are equivalent to shifting covariates.

    The covariance hyperparameter can be specified either as a vector of standard
    deviations, using the \code{sd} argument, a vector of variances using the
    \code{cov} argument, or the entire variance/covariance matrix itself. When
    specifying standard deviations, a vector of length less than the number of fixed effects will
    have its tail repeated, while the first element is assumed to apply only
    to the intercept term. So in the default of \code{c(10, 2.5)}, the intercept
    receives a standard deviation of 10 and the various slopes are all given
    a standard deviation of 2.5.

    The \code{common.scale} argument specifies whether or not the
    prior is to be interpretted as being on the same scale as the residuals.
    To specify a prior in an absolute sense, set to \code{FALSE}. Argument
    is only applicable to linear mixed models.

    \item \code{t(df = 3, mean = 0, scale = c(10^2, 2.5^2), common.scale = TRUE)}
    
    The degrees of freedom - \code{df} argument - must be positive. If \code{mean} is
    of length 1, it is repeated for every fixed effect. Length 2 repeats just the second
    element for all slopes. Otherwise, the length must be equal to that of the number of
    fixed effects.
    
    If \code{scale} is of length 1, it is repeated along the diagonal for every
    component. Length 2 repeats just the second element for all slopes. Length equal
    to the number of fixed effects sees the vector simply turned into a diagonal matrix.
    Finally, it can be a full scale matrix, so long as it is positive definite.

    \code{t} priors for linear mixed models require that the fixed effects be added to
    set of parameters that are numerically optimized, and thus can substantially
    increase running time. In addition, when \code{common.scale} is \code{TRUE}, the
    residual variance must be numerically optimized as well. \code{normal} priors
    on the common scale can be fully profiled and do not suffer from this drawback.

    At present, \code{t} priors cannot be used with the \code{REML = TRUE} argument
    as that implies an integral without a closed form solution.

    \item \code{horseshoe(mean = 0, global.shrinkage = 2.5, common.scale = TRUE)}
    
    The horseshoe shrinkage prior is implemented similarly to the \code{t} prior, in
    that it requires adding the fixed effects to the parameter set for numeric optimization.
    
    \code{global.shrinkage}, also referred to as \eqn{\tau}, must be positive and
    is on the scale of a standard deviation. Local shrinkage parameters are treated as
    independent across all fixed effects and are integrated out. See
    \cite{Carvalho et al. (2009)} in the references.
  }

  \strong{Covariance Priors}

  \itemize{
    \item \code{gamma(shape = 2.5, rate = 0, common.scale = TRUE, posterior.scale = "sd")}

    Applicable only for univariate grouping factors. A
    rate of \code{0} or a shape of \code{0} imposes an improper prior. The
    posterior scale can be \code{"sd"} or \code{"var"} and determines the scale
    on which the prior is meant to be applied.

    \item \code{invgamma(shape = 0.5, scale = 10^2, common.scale = TRUE, posterior.scale = "sd")}
    
    Applicable only for univariate grouping factors. A
    scale of \code{0} or a shape of \code{0} imposes an improper prior. Options
    are as above.

   \item \code{wishart(df = level.dim + 2.5, scale = Inf, common.scale = TRUE, posterior.scale = "cov")}

   A scale of \code{Inf} or a shape of \code{0} imposes an improper prior. The behavior
   for singular matrices with only some infinite eigenvalues is undefined. Posterior scale
   can be \code{"cov"} or \code{"sqrt"},
   the latter of which applies to the unique matrix root that is also a valid covariance
   matrix.

   \item \code{invwishart(df = level.dim - 0.5, scale = diag(10^2 / (df + level.dim + 1), level.dim),
                        common.scale = TRUE, posterior.scale = "cov")}

   A scale of \code{0} or a shape of \code{0} imposes an improper prior. The behavior
   for singular matrices with only some zero eigenvalues is undefined.

   \item \code{custom(fn, chol = FALSE, common.scale = TRUE, scale = "none")}

   Applies to the given function (\code{fn}). If \code{chol} is \code{TRUE}, \code{fn} is
   passed a \emph{right} factor of covariance matrix; \code{FALSE} results in the matrix being
   passed directly. \code{scale} can be \code{"none"}, \code{"log"}, or \code{"dev"} corresponding to
   \eqn{p(\Sigma)}, \eqn{\log p(\Sigma)}{log p(\Sigma)}, and \eqn{-2 \log p(\Sigma)}{-2 log p(\Sigma)}
   respectively.

   Since the prior is may have an arbitrary form, setting \code{common.scale} to \code{FALSE}
   for a linear mixed model means that full profiling may no longer be possible. As such,
   that parameter is numerically optimized.
  }

  \strong{Residual Variance Priors}

  \itemize{
    \item \code{point(value = 1.0, posterior.scale = "sd")}
  
    Fixes the parameter to a specific value given as either an \code{"sd"} or a \code{"var"}.

    \item \code{gamma(shape = 0, rate = 0, posterior.scale = "var")}
    
    As above with different defaults.

    \item \code{invgamma(shape = 0, scale = 0, posterior.scale = "var")}

    As above with different defaults.
  }
}

\section{Evaluating Environment}{
  The variables that the defining environment have populated are:
  \itemize{
    \item \code{p} aliased to \code{n.fixef} - the number of fixed effects
    \item \code{n} aliased to \code{n.obs} - the number of observations
    \item \code{q.k} aliased to \code{level.dim} - for covariance priors, the dimension of the grouping factor/grouping level
    \item \code{j.k} aliased to \code{n.grps} - also for covariance priors, the number of groups that comprise a specific grouping factor
  }
}

\section{Methods}{
  \describe{
    \item{toString}{Pretty-prints the distribution and its parameters.}
  }
}

\references{
Carvalho, Carlos M., Nicholas G. Polson, and James G. Scott.
  "Handling Sparsity via the Horseshoe."
  AISTATS. Vol. 5. 2009.
}

\seealso{
  \code{\link{blmer}()} and \code{\link{bglmer}()},
  which produce these objects, and \code{\link{bmerMod-class}} objects which contain them.
}
\keyword{classes}
