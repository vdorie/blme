setClass("bmerHorseshoeDist",
         representation(beta.0 = "numeric",
                        tau.sq = "numeric"),
         contains = "bmerDist")

toString.bmerHorseshoeDist <- function(x, digits = getOption("digits"), ...) {  
  meanString <- ""
  beta.0 <- x@beta.0
  if (length(beta.0) > 4L) {
    meanString <- paste0("mean = c(", toString(round(beta.0[seq_len(4L)], digits)), ", ...)")
  } else if (length(beta.0) == 1L) {
    meanString <- paste0("mean = ", toString(round(beta.0[1L], digits)))
  } else {
    meanString <- paste0("mean = c(", toString(round(beta.0, digits)), ")")
  }
     
  paste0("horseshoe(", meanString, ", ",
         "global.shrinkage = ", round(sqrt(x@tau.sq), digits), ", ",
         "common.scale = ", x@commonScale, ")")
}
setMethod("getDFAdjustment", "bmerHorseshoeDist",
  function(object) {
    if (object@commonScale == TRUE) length(object@beta.0) else 0
  }
)
setMethod("getConstantTerm", "bmerHorseshoeDist",
  function(object) {
    d <- length(object@beta.0)
    
    d * (3 * log(pi) + log(2) + log(object@tau.sq))
  }
)
setMethod("getExponentialTerm", "bmerHorseshoeDist",
  function(object, beta, sigma = NULL) {
    beta.0 <- object@beta.0
    tau.sq <- object@tau.sq
    
    dist <- 0.5 * (beta - beta.0)^2 / tau.sq
    if (object@commonScale == TRUE && !is.null(sigma)) dist <- dist / sigma^2
    
    temp <- suppressWarnings(sapply(dist, expint::expint_E1, scale = TRUE))
    temp[is.nan(temp)] <- .Machine$double.xmax * .Machine$double.eps
    
    result <- -2 * sum(log(temp))
    
    c(0, result)
  }
)

