setClass("bmerTDist",
         representation(df = "numeric",
                        R.scale.inv = "matrix"),
         contains = "bmerDist")

toString.bmerTDist <- function(x, digits = getOption("digits"), ...) {
  scaleString <- ""
  scale <- crossprod(solve(x@R.scale.inv))
  
  if (nrow(scale) > 2) {
    scaleString <- paste("scale = c(", toString(round(scale[1:4], digits)), ", ...)", sep = "")
  } else if (nrow(scale) == 2) {
    scaleString <- paste("scale = c(", toString(round(scale[1:4], digits)), ")", sep = "")
  } else {
    scaleString <- paste("scale = ", toString(round(scale[1], digits)), sep = "")
  }
  
  paste("t(df = ", x@df, ", ", scaleString,
        ", common.scale = ", x@commonScale,
        ")", sep="")
}
setMethod("getDFAdjustment", "bmerTDist",
  function(object) {
    if (object@commonScale == TRUE) nrow(object@R.scale.inv) else 0
  }
)
setMethod("getConstantTerm", "bmerTDist",
  function(object) {
    R.scale.inv <- object@R.scale.inv
    d <- nrow(R.scale.inv)
    df <- object@df
    
    -2.0 * lgamma(0.5 * (df + d)) + 2.0 * lgamma(0.5 * df) +
      d * (log(df) + log(pi)) - 2.0 * sum(log(diag(R.scale.inv)))
  }
)
setMethod("getExponentialTerm", "bmerTDist",
  function(object, beta) {
    R.scale.inv <- object@R.scale.inv
    d <- nrow(R.scale.inv)
    df <- object@df

    dist <- tcrossprod(crossprod(beta, R.scale.inv))[1]
    
    exponential <- (df + d) * log(1 + dist / df)
    c(0, exponential)
  }
)

