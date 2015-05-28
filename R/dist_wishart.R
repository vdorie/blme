setClass("bmerWishartDist",
         representation(df          = "numeric",
                        R.scale.inv = "matrix",
                        log.det.scale = "numeric",
                        posteriorScale = "character"),
         contains = "bmerDist",
         validity = function(object) object@posteriorScale == "cov" || object@posteriorScale == "sqrt")

setClass("bmerInvWishartDist",
         representation(df      = "numeric",
                        R.scale = "matrix",
                        log.det.scale = "numeric",
                        posteriorScale = "character"),
         contains = "bmerDist",
         validity = function(object) object@posteriorScale == "cov" || object@posteriorScale == "sqrt")

toString.bmerWishartDist <- function(x, digits = getOption("digits"), ...) {
  if (any(diag(x@R.scale.inv) == 0)) {
    scale <- Inf
  } else if (any(is.infinite(x@R.scale.inv))) {
    scale <- 0
  } else {
    scale <- solve(tcrossprod(x@R.scale.inv))
  }
  
  if (length(scale) == 1) {
    scaleString <- round(scale, digits)
  } else if (nrow(scale) > 2) {
    scaleString <- paste("c(", toString(round(scale[1:3], digits)), ", ...)", sep = "")
  } else if (nrow(scale) == 2) {
    scaleString <- paste("c(", toString(round(scale[1:4], digits)), ")", sep = "")
  }
  
  paste("wishart(",
        "df = ", round(x@df, digits),
        ", scale = ", scaleString,
        ", posterior.scale = ", x@posteriorScale,
        ", common.scale = ", x@commonScale,
        ")", sep="")
}

toString.bmerInvWishartDist <- function(x, digits = getOption("digits"), ...) {
  if (any(diag(x@R.scale) == 0)) {
    scale <- 0
  } else if (any(is.infinite(x@R.scale))) {
    scale <- Inf
  } else {
    scale <- crossprod(x@R.scale)
  }
  
  if (length(scale) == 1) {
    scaleString <- round(scale, digits)
  } else if (nrow(scale) > 2) {
    scaleString <- paste("c(", toString(round(scale[1:3], digits)), ", ...)", sep = "")
  } else if (nrow(scale) == 2) {
    scaleString <- paste("c(", toString(round(scale[1:4], digits)), ")", sep = "")
  }
  
  paste("invwishart(",
        "df = ", round(x@df, digits),
        ", scale = ", scaleString,
        ", posterior.scale = ", x@posteriorScale,
        ", common.scale = ", x@commonScale,
        ")", sep="")
}

setMethod("getDFAdjustment", "bmerWishartDist",
  function(object) {
    factorDim <- nrow(object@R.scale.inv)
    if (object@commonScale || !is.finite(object@log.det.scale)) 0 else -factorDim * (object@df - factorDim - 1.0)
  }
)

setMethod("getDFAdjustment", "bmerInvWishartDist",
  function(object) {
    factorDim <- nrow(object@R.scale)
    if (object@commonScale || !is.finite(object@log.det.scale)) 0 else factorDim * (object@df + factorDim + 1.0)
  }
)

setMethod("getConstantTerm", "bmerWishartDist",
  function(object) {
    df <- object@df; R.scale.inv <- object@R.scale.inv
    log.det.scale <- object@log.det.scale
    
    if (is.infinite(log.det.scale)) return (0.0)

    factorDim <- nrow(R.scale.inv)
    
    result <- df * (factorDim * log(2) + log.det.scale) +
      0.5 * factorDim * (factorDim - 1.0) * log(pi)
    for (i in 1:factorDim)
      result <- result + 2.0 * lgamma(0.5 * (df + 1.0 - i))

    result
  }
)

setMethod("getConstantTerm", "bmerInvWishartDist",
  function(object) {
    df <- object@df; R.scale <- object@R.scale
    log.det.scale <- object@log.det.scale
    
    if (is.infinite(log.det.scale)) return (0.0)

    factorDim <- nrow(R.scale)

    result <- df * (factorDim * log(2) - log.det.scale) +
      0.5 * factorDim * (factorDim - 1.0) * log(pi)
    for (i in 1:factorDim)
      result <-  result + 2.0 * lgamma(0.5 * (df + 1.0 - i))

    result
  }
)

setMethod("getExponentialSigmaPower", "bmerWishartDist",
  function (object) {
    if (object@commonScale) return(0)
    
    if (object@posteriorScale == "sqrt") 1 else 2
  })

setMethod("getExponentialSigmaPower", "bmerInvWishartDist",
  function (object) {
    if (object@commonScale) return(0)
    
    if (object@posteriorScale == "sqrt") -1 else -2
  })


setMethod("getExponentialTerm", "bmerWishartDist",
  function(object, Lambda.t) {
    if (is.infinite(object@log.det.scale)) return(c(0, 0.0))

    if (object@posteriorScale == "cov") {
      temp <- Lambda.t %*% object@R.scale.inv
      exponential <- sum(temp^2)
      
      power <- 2
    } else {
      Sigma <- crossprod(Lambda.t)
      decomp <- eigen(Sigma)
      Sigma.sqrt <- decomp$vectors %*% tcrossprod(diag(sqrt(decomp$values)), decomp$vectors)
      exponential <- sum(Sigma.sqrt * crossprod(object@R.scale.inv))

      power <- 1
    }
    
    if (object@commonScale) c(0, exponential) else c(power, exponential)
  }
)

setMethod("getExponentialTerm", "bmerInvWishartDist",
  function(object, Lambda.t) {
    if (is.infinite(object@log.det.scale)) return(c(0, 0.0))
    
    if (object@posteriorScale == "cov") {
      power <- -2
      
      if (any(diag(Lambda.t) == 0))
        return (if (object@commonScale) c(0, Inf) else c(power, Inf))
        
      temp <- object@R.scale %*% solve(Lambda.t)
      exponential <- sum(temp^2)

    } else {
      power <- -1
      
      if (any(diag(Lambda.t) == 0))
        return (if (object@commonScale) c(0, Inf) else c(power, Inf))
      
      Sigma <- crossprod(Lambda.t)
      decomp <- eigen(Sigma)
      Sigma.inv.sqrt <- decomp$vectors %*% tcrossprod(diag(1 / sqrt(decomp$values)), decomp$vectors)
      exponential <- sum(Sigma.inv.sqrt * tcrossprod(object@R.scale))
    }

    if (object@commonScale) c(0, exponential) else c(power, exponential)
  }
)

setMethod("getPolynomialTerm", "bmerWishartDist",
  function(object, Lambda.t) {
    factorDim <- nrow(object@R.scale.inv)
    -2.0 * (object@df - factorDim - 1.0) * sum(log(diag(Lambda.t)))
  }
)

setMethod("getPolynomialTerm", "bmerInvWishartDist",
  function(object, Lambda.t) {
    factorDim <- nrow(object@R.scale)
   2.0 * (object@df + factorDim + 1.0) * sum(log(diag(Lambda.t)))
  }
)
