\name{bmer-class}
\title{Bayesian Linear Mixed-Effects Model Representations and *bmer Methods}
\docType{class}
% Classes
\alias{bmer-class}
\alias{summary.bmer-class}
% Methods
\alias{print,bmer-method}
\alias{show,bmer-method}
\alias{summary,bmer-method}

\description{The \code{bmer} class represents linear or generalized
  linear or nonlinear mixed-effects models with possible priors over
  model components. It inherits from the \code{\linkS4class{mer}}
  class.
}

\section{Objects from the Class}{Objects can be created by calls of the
  form \code{new("bmer", ...)}  or, more commonly, via the
  \code{\link{blmer}} or \code{\link{bglmer}} functions.
}

\section{Slots}{
  The class \code{"bmer"} contains three slots beyond that of the
  \code{"mer"} class:

  \describe{
    \item{\code{cov.prior}:}{Prior over the covariance of the 
      modeled coefficients, or "random effects".}
    \item{\code{fixef.prior}:}{Prior over the unmodeled coefficients,
      or "fixed effects".}
    \item{\code{var.prior}:}{Prior over the observation level
      variance, or "common scale" parameter.}
    \item{\code{A}, \code{Cm}, \code{Cx}, \code{Gp}, \code{L},
          \code{RX}, \code{RZX}, \code{ST}, \code{V}, \code{X},
          \code{Xst}, \code{Zt}, \code{approximate}, \code{call},
          \code{deviance}, \code{dims}, \code{env}, \code{eta},
          \code{fixef}, \code{flist}, \code{frame}, \code{ghw},
          \code{ghx}, \code{mu}, \code{muEta}, \code{nlmodel},
          \code{offset}, \code{pWt}, \code{ranef}, \code{resid},
          \code{sqrtXWt}, \code{sqrtrWt}, \code{u}, \code{var},
          \code{y}}{See \code{\linkS4class{mer}}.}
  }
}
\seealso{
  \code{\link{blmer}()} and \code{\link{bglmer}()},
  which produce these objects.
  \cr
  \code{\linkS4class{mer}}, from which this class inherits.
}
\keyword{classes}
