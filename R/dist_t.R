setClass("bmerTDist",
         representation(df = "numeric",
                        beta.0 = "numeric",
                        R.scale.inv = "matrix"),
         contains = "bmerDist")

toString.bmerTDist <- function(x, digits = getOption("digits"), ...) {
  scaleString <- ""
  scale <- crossprod(solve(x@R.scale.inv))
    
  meanString <- ""
  beta.0 <- x@beta.0
  if (length(beta.0) > 4) {
    meanString <- paste0("mean = c(", toString(round(beta.0[seq_len(4)], digits)), ", ...)")
  } else if (length(beta.0) == 1) {
    meanString <- paste0("mean = ", toString(round(beta.0[1], digits)))
  } else {
    meanString <- paste0("mean = c(", toString(round(beta.0[seq_len(4)], digits)), ")")
  }
  
  if (nrow(scale) > 2) {
    scaleString <- paste0("scale = c(", toString(round(scale[seq_len(4)], digits)), ", ...)")
  } else if (nrow(scale) == 2) {
    scaleString <- paste0("scale = c(", toString(round(scale[seq_len(4)], digits)), ")")
  } else {
    scaleString <- paste0("scale = ", toString(round(scale[1], digits)))
  }
    
  paste("t(df = ", x@df, ", ", meanString, ", ", scaleString,
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
    beta.0 <- object@beta.0
    R.scale.inv <- object@R.scale.inv
    d <- nrow(R.scale.inv)
    df <- object@df

    dist <- tcrossprod(crossprod(beta - beta.0, R.scale.inv))[1]
    
    exponential <- (df + d) * log(1 + dist / df)
    c(0, exponential)
  }
)

