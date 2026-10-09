setClass("bmerIndependentTsDist",
         representation(df = "numeric",
                        beta.0 = "numeric",
                        d = "numeric",
                        sqrt.scale.inv = "numeric"),
         contains = "bmerDist")

toString.bmerIndependentTsDist <- function(x, digits = getOption("digits"), ...) {
  vectorToString <- function(name, values) {
    if (length(values) > 1L && all(values == values[1L])) values <- values[1L]
    if (length(values) > 4L) {
      paste0(name, " = c(", toString(round(values[seq_len(4L)], digits)), ", ...)")
    } else if (length(values) == 1L) {
      paste0(name, " = ", toString(round(values[1L], digits)))
    } else {
      paste0(name, " = c(", toString(round(values, digits)), ")")
    }
  }

  dfString    <- vectorToString("df", x@df)
  meanString  <- vectorToString("mean", x@beta.0)
  scaleString <- vectorToString("scale", 1 / x@sqrt.scale.inv^2)

  paste0("independent_ts(", dfString, ", ", meanString, ", ", scaleString, ", ",
         "common.scale = ", x@commonScale, ")")
}
setMethod("getDFAdjustment", "bmerIndependentTsDist",
  function(object) {
    if (object@commonScale == TRUE) object@d else 0
  }
)
setMethod("getConstantTerm", "bmerIndependentTsDist",
  function(object) {
    ## components with infinite scale are flat and contribute nothing
    finite <- object@sqrt.scale.inv > 0
    sqrt.scale.inv <- object@sqrt.scale.inv[finite]
    df <- object@df[finite]

    -2.0 * sum(lgamma(0.5 * (df + 1))) + 2.0 * sum(lgamma(0.5 * df)) +
      sum(log(df) + log(pi)) - 2.0 * sum(log(sqrt.scale.inv))
  }
)
setMethod("getExponentialTerm", "bmerIndependentTsDist",
  function(object, beta, sigma = NULL) {
    finite <- object@sqrt.scale.inv > 0
    beta.0 <- object@beta.0[finite]
    sqrt.scale.inv <- object@sqrt.scale.inv[finite]
    df <- object@df[finite]

    dist <- ((beta[finite] - beta.0) * sqrt.scale.inv)^2
    if (object@commonScale == TRUE && !is.null(sigma)) {
      ## density is 0 at sigma = 0; sigma is bounded away from it during optimization
      if (sigma <= 0) return(c(0, Inf))
      dist <- dist / sigma^2
    }
    if (any(is.na(dist)) || any(is.infinite(dist))) stop("non-finite or NA result in independent_ts prior")

    exponential <- sum((df + 1) * log(1 + dist / df))
    c(0, exponential)
  }
)
