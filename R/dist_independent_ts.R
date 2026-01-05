setClass("bmerIndependentTsDist",
         representation(df = "numeric",
                        beta.0 = "numeric",
                        d = "numeric",
                        sqrt.scale.inv = "numeric"),
         contains = "bmerDist")

toString.bmerTDist <- function(x, digits = getOption("digits"), ...) {
  dfString <- ""
  if (length(df) > 4L) {
    meanString <- paste0("df = c(", toString(round(df[seq_len(4L)], digits)), ", ...)")
  } else if (length(beta.0) == 1L) {
    meanString <- paste0("df = ", toString(round(df[1L], digits)))
  } else {
    meanString <- paste0("df = c(", toString(round(df, digits)), ")")
  }

  meanString <- ""
  beta.0 <- x@beta.0
  if (length(beta.0) > 4L) {
    meanString <- paste0("mean = c(", toString(round(beta.0[seq_len(4L)], digits)), ", ...)")
  } else if (length(beta.0) == 1L) {
    meanString <- paste0("mean = ", toString(round(beta.0[1L], digits)))
  } else {
    meanString <- paste0("mean = c(", toString(round(beta.0, digits)), ")")
  }

  scaleString <- ""
  scale <- 1 / scale.inv^2
  if (length(scale) > 4L) {
    scaleString <- paste0("scale = c(", toString(round(scale[seq_len(4L)], digits)), ", ...")
  } else if (length(scale) == 1L) {
    scaleString <- paste0("scale = ", toString(round(scale[1L], digits)))
  } else {
    scaleString <- paste0("scale = c(", toString(round(scale, digits)), ")")
  }
    
  paste0("independent_ts(df = ", dfString, ", ", meanString, ", ", scaleString, ", ",
         "common.scale = ", x@commonScale, ")")
}
setMethod("getDFAdjustment", "bmerIndependentTsDist",
  function(object) {
    if (object@commonScale == TRUE) object@d else 0
  }
)
setMethod("getConstantTerm", "bmerIndependentTsDist",
  function(object) {
    sqrt.scale.inv <- object@sqrt.scale.inv
    d <- object@d
    df <- object@df
    

    ldet <- sum(log(sqrt.scale.inv[sqrt.scale.inv > 0]))
    
    -2.0 * sum(lgamma(0.5 * (df + 1))) + 2.0 * sum(lgamma(0.5 * df)) +
      sum(log(df) + log(pi)) - 2.0 * ldet
  }
)
setMethod("getExponentialTerm", "bmerIndependentTsDist",
  function(object, beta, sigma = NULL) {
    beta.0 <- object@beta.0
    sqrt.scale.inv <- object@sqrt.scale.inv
    d <- object@d
    df <- object@df

    dist <- ((beta - beta.0) * sqrt.scale.inv)^2
    if (object@commonScale == TRUE && !is.null(sigma)) dist <- dist / sigma^2
    if (any(is.na(dist)) || any(is.infinite(dist))) stop("non-finite or NA result in independent ts-prior")
    
    exponential <- sum((df + 1) * log(1 + dist / df))
    c(0, exponential)
  }
)

