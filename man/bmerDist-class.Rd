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
  form \code{new("bmerPrior", ...)}  or, more commonly, via the
  \code{\link{blmer}} or \code{\link{bglmer}} functions.

  When using the main \code{blme} functions, the prior-related arguments can be
  passed what essentially
  are function calls with the distinction that they are delayed in evaluation
  until information about the model is available. At that time, the functions
  are defined in a special environment and then evaluated in the calling
  environment. This is reflected in the prototypes of various prior-creating
  functions which depend on parameters not available in the top-level
  environment.

  Finally, if the trailing parentheses are omitted from a \code{blmer}/\code{bglmer}
  prior argument, they are simply added as a form of \dQuote{syntactic sugar}.
}

\section{Prior Distributions}{
  This section lists the prototypes for the functions that are called to parse
  a prior during a model fit.

  \strong{Fixed Effect Priors}

  \itemize{
    \item \code{normal(sd = c(10, 2.5), cov, common.scale = TRUE)}
    
    Applies to \code{fixef.prior} argument. Gaussian prior for the fixed
    effects are currently constrained to have a mean of 0. The covariance
    can be specified either as a vector of standard deviations, using the 
    \code{sd} argument, a vector of variances using the \code{cov} argument,
    or the entire variance/covariance matrix itself. When specifying standard
    deviations, a vector of length less than the number of fixed effects will
    have its tail repeated, while the first element is assumed to apply only
    to the intercept term. So in the default of \code{c(10, 2.5)}, the intercept
    receives a standard deviation of 10 and the various slopes are all given
    a standard deviation of 2.5.

    The \code{common.scale} argument specifies whether or not the
    prior is to be interpretted as being on the same scale as the residuals.
    To specify a prior in an absolute sense, set to \code{FALSE}. Argument
    is only present for linear mixed models.
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
  The variables that the evaluating environment have populated are:
  \itemize{
    \item \code{p} aliased to \code{n.fixef} - the number of fixed effects
    \item \code{n} aliased to \code{n.obs} - the number of observations
    \item \code{level.dim} - for covariance priors, the dimension of the grouping factor/grouping level
  }
}

\section{Methods}{
  \describe{
    \item{toString}{Pretty-prints the distribution and its parameters.}
  }
}

\seealso{
  \code{\link{blmer}()} and \code{\link{bglmer}()},
  which produce these objects, and \code{\link{bmerMod-class}} objects which contain them.
}
\keyword{classes}
