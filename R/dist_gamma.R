setClass("bmerGammaDist",
         representation(shape = "numeric",
                        rate  = "numeric",
                        posteriorScale = "character"),
         contains = "bmerDist",
         validity = function(object) object@posteriorScale == "var" || object@posteriorScale == "sd")
setClass("bmerInvGammaDist",
         representation(shape = "numeric",
                        scale = "numeric",
                        posteriorScale = "character"),
         contains = "bmerDist",
         validity = function(object) object@posteriorScale == "var" || object@posteriorScale == "sd")

toString.bmerGammaDist <- function(x, digits = getOption("digits"), includeCommonScale = TRUE, ...) {
  paste("gamma(shape = ", round(x@shape, digits),
        ", rate = ", round(x@rate, digits),
        ", posterior.scale = ", x@posteriorScale,
        if (includeCommonScale) paste(", common.scale = ", x@commonScale, sep = "") else "",
        ")", sep="")
}

toString.bmerInvGammaDist <- function(x, digits = getOption("digits"), includeCommonScale = TRUE, ...) {
  paste("invgamma(shape = ", round(x@shape, digits),
        ", scale = ", round(x@scale, digits),
        ", posterior.scale = ", x@posteriorScale,
        if (includeCommonScale) paste(", common.scale = ", x@commonScale, sep = "") else "",
        ")", sep="")
}

setMethod("getDFAdjustment", "bmerGammaDist",
  function(object) {
    if (object@commonScale) 0 else {
      if (object@posteriorScale == 'sd')
             -(object@shape - 1.0)
      else
        -2.0 * (object@shape - 1.0)
    }
  }
)
setMethod("getDFAdjustment", "bmerInvGammaDist",
  function(object) {
    if (object@commonScale) 0 else {
      if (object@posteriorScale == 'sd')
              (object@shape + 1.0)
      else
        2.0 * (object@shape + 1.0)
    }
  }
)

setMethod("getConstantTerm", "bmerGammaDist",
  function(object) {
    shape <- object@shape; rate <- object@rate
    if (shape == 0.0 || rate == 0.0) return(0.0)
    if (shape <  0.0 || rate <  0.0) return(NaN)

    -2.0 * (shape * log(rate) - lgamma(shape))
  }
)
setMethod("getConstantTerm" ,"bmerInvGammaDist",
  function(object) {
    shape <- object@shape; scale <- object@scale
    
    if (shape == 0.0 || scale == 0.0) return(0.0)
    if (shape <  0.0 || scale <  0.0) return(NaN)
  
    -2.0 * (shape * log(scale) - lgamma(shape))
  }
)

setMethod("getExponentialSigmaPower", "bmerGammaDist",
  function (object) {
    if (object@commonScale || object@rate == 0) return(0)
    
    if (object@posteriorScale == "sd") 1 else 2
  })

setMethod("getExponentialSigmaPower", "bmerInvGammaDist",
  function (object) {
    if (object@commonScale || object@scale == 0) return(0)
    
    if (object@posteriorScale == "sd") -1 else -2
  })

setMethod("getExponentialTerm", "bmerGammaDist",
  function(object, lambda) {
    if (object@rate == 0) return (c(0, 0.0))

    if (missing(lambda)) lambda <- 1.0
    
    if (object@posteriorScale == "var") {
      exponential <- 2.0 * lambda^2 * object@rate
      power <- 2
    } else {
      exponential <- 2.0 * lambda * object@rate
      power <- 1
    }

    if (object@commonScale == TRUE) c(0, exponential) else c(power, exponential)
  })

setMethod("getExponentialTerm", "bmerInvGammaDist",
  function(object, lambda) {
    if (object@scale == 0) return (c(0, 0.0))

    if (missing(lambda)) lambda <- 1.0
    
    if (object@posteriorScale == "var") {
      exponential <- 2.0 * object@scale / lambda^2
      power <- -2
    } else {
      exponential <- 2.0 * object@scale / lambda
      power <- -1
    }

    if (object@commonScale == TRUE) c(0, exponential) else c(power, exponential)
  })

setMethod("getPolynomialTerm", "bmerGammaDist",
  function(object, lambda) {
    if (object@posteriorScale == "var")
      -4 * (object@shape - 1.0) * log(lambda)
    else
      -2 * (object@shape - 1.0) * log(lambda)
  }
)

setMethod("getPolynomialTerm", "bmerInvGammaDist",
  function(object, lambda) {
    if (object@posteriorScale == "var")
      4.0 * (object@shape + 1.0) * log(lambda)
    else
      2.0 * (object@shape + 1.0) * log(lambda)
  }
)
