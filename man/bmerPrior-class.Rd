\name{bmerPrior-class}
\title{Bayesian Linear Mixed-Effects Model Prior Representations and *bmerPrior Methods}
\docType{class}
% Classes
\alias{bmerPrior-class}
% Methods
\alias{print,bmerPrior-method}
\alias{show,bmerPrior-method}

\description{Prints a \code{bmerPrior} object in a human-readable format. Mostly used in debugging.}

\section{Objects from the Class}{Objects can be created by calls of the
  form \code{new("bmerPrior", ...)}  or, more commonly, via the
  \code{\link{blmer}} or \code{\link{bglmer}} functions.
}

\section{Slots}{
  The class \code{"bmerPrior"} contains four slots:

  \describe{
    \item{\code{type}:}{An integer representing which "type" of prior it is,
          		taken from an enumeration. Types refer to how the
			parameter space is to be decomposed, if at all.}
    \item{\code{families}:}{An integer-enum of which parametric families the
                            prior contains. For some decompositions, this is
                            a vector. For most, it is a single value.}
    \item{\code{scales}:}{An integer-enum detailing the scale of the desired
                          posterior, typically chosen between standard
                          deviation and variance.}
    \item{\code{hyperparameters}:}{How the family is parameterized. Values
          are typically stored in a way that facilitates computation, and may
          not be easily interpretted.}
  }
}
\seealso{
  \code{\link{blmer}()} and \code{\link{bglmer}()},
  which produce these objects, and \code{\link{bmer-class}} objects which contain them.
}
\keyword{classes}
