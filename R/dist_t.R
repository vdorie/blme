setClass("bmerTDist",
         representation(df = "numeric",
                        beta.0 = "numeric",
                        d = "numeric",
                        R.scale.inv = "matrix"),
         contains = "bmerDist")

toString.bmerTDist <- function(x, digits = getOption("digits"), ...) {
  scaleString <- ""
  scale <- if (x@d != nrow(x@R.scale.inv)) diag(1 / diag(x@R.scale.inv)^2) else crossprod(solve(x@R.scale.inv))
    
  meanString <- ""
  beta.0 <- x@beta.0
  if (length(beta.0) > 4L) {
    meanString <- paste0("mean = c(", toString(round(beta.0[seq_len(4L)], digits)), ", ...)")
  } else if (length(beta.0) == 1L) {
    meanString <- paste0("mean = ", toString(round(beta.0[1L], digits)))
  } else {
    meanString <- paste0("mean = c(", toString(round(beta.0, digits)), ")")
  }
  if (all(scale[lower.tri(scale) | upper.tri(scale)] == 0)) {
    scale <- diag(scale)
    if (length(scale) > 4L) {
      scaleString <- paste0("diag(scale) = c(", toString(round(scale[seq_len(4L)], digits)), ", ...")
    } else if (length(scale) == 1L) {
      scaleString <- paste0("scale = ", toString(round(scale[1L], digits)))
    } else {
      scaleString <- paste0("diag(scale) = c(", toString(round(scale, digits)), ")")
    }
  } else {
    if (nrow(scale) > 2) {
      scaleString <- paste0("scale = c(", toString(round(scale[seq_len(4)], digits)), ", ...)")
    } else if (nrow(scale) == 2) {
      scaleString <- paste0("scale = c(", toString(round(scale[seq_len(4)], digits)), ")")
    } else {
      scaleString <- paste0("scale = ", toString(round(scale[1], digits)))
    }
  }
    
  paste("t(df = ", x@df, ", ", meanString, ", ", scaleString,
        ", common.scale = ", x@commonScale,
        ")", sep="")
}
setMethod("getDFAdjustment", "bmerTDist",
  function(object) {
    if (object@commonScale == TRUE) object@d else 0
  }
)
setMethod("getConstantTerm", "bmerTDist",
  function(object) {
    R.scale.inv <- object@R.scale.inv
    d <- object@d
    df <- object@df
    
    det <- sum(log(if (d != nrow(R.scale.inv)) { p <- diag(R.scale.inv); p[p > 0] } else diag(R.scale.inv)))
    
    -2.0 * lgamma(0.5 * (df + d)) + 2.0 * lgamma(0.5 * df) +
      d * (log(df) + log(pi)) - 2.0 * det
  }
)
setMethod("getExponentialTerm", "bmerTDist",
  function(object, beta) {
    beta.0 <- object@beta.0
    R.scale.inv <- object@R.scale.inv
    d <- object@d
    df <- object@df

    dist <- tcrossprod(crossprod(beta - beta.0, R.scale.inv))[1L]
    if (any(is.na(dist)) || any(is.infinite(dist))) browser()
    
    exponential <- (df + d) * log(1 + dist / df)
    c(0, exponential)
  }
)

