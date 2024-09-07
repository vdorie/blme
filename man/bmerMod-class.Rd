\docType{class}
\name{bmerMod-class}
\title{Class "bmerMod" of Fitted Mixed-Effect Models} 
\alias{blmerMod-class}
\alias{bglmerMod-class}
\alias{bmerMod}
\alias{bmerMod-class}
\alias{print.bmerMod}
\alias{print.summary.bmerMod}
\alias{show,bmerMod-method}
\alias{show.bmerMod}
\alias{summary.bmerMod}
\alias{summary.summary.bmerMod}
\alias{vcov.summary.bmerMod}

\description{The \code{bmerMod} class represents linear or generalized
  linear or nonlinear mixed-effects models with possible priors over
  model components. It inherits from the \code{\link[lme4:merMod-class]{merMod}}
  class.
}

\section{Objects from the Class}{
  Objects are created by calls to \code{\link{blmer}} or \code{\link{bglmer}}.
}

\section{Slots}{
  A \code{bmerMod} object contains one additional slot beyond the base
  \code{merMod} class:

  \describe{
    \item{\code{priors}:}{A named list comprised of \code{covPriors},
      \code{fixefPrior}, and \code{residPrior}.}
  }
 
  In addition, the \code{devcomp} slot, element \code{cmp} includes the
  \code{penalty} item which is the computed deviance for the priors. Add
  this to the regular deviance to obtain the value of the objective function
  that is used in optimization.
}
\seealso{
  \code{\link{blmer}} and \code{\link{bglmer}},
  which produce these objects.
  \cr
  \code{\link[lme4:merMod-class]{merMod}}, from which this class inherits.
}
\examples{
showClass("bmerMod")
methods(class = "bmerMod")
}
\keyword{classes}
