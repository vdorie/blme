setClass("bmerDist", representation(commonScale = "logical"),
         contains = "VIRTUAL")

if (!isGeneric("getDFAdjustment"))
  setGeneric("getDFAdjustment", function(object, ...) standardGeneric("getDFAdjustment"))
if (!isGeneric("getConstantTerm"))
  setGeneric("getConstantTerm", function(object, ...) standardGeneric("getConstantTerm"))
## what power sigma has if prior induces exp(-0.5 * sigma^pow * stuff)
if (!isGeneric("getExponentialSigmaPower"))
  setGeneric("getExponentialSigmaPower", function(object, ...) standardGeneric("getExponentialSigmaPower"))
## whatever is going on in the exponent, and what power sigma has connected with it
## note, relative to 1 / 2
if (!isGeneric("getExponentialTerm"))
  setGeneric("getExponentialTerm", function(object, ...) standardGeneric("getExponentialTerm"))
if (!isGeneric("getPolynomialTerm"))
  setGeneric("getPolynomialTerm", function(object, ...) standardGeneric("getPolynomialTerm"))

setMethod("getDFAdjustment", "ANY", function(object, ...) 0)
setMethod("getConstantTerm", "ANY", function(object, ...) 0)
setMethod("getExponentialSigmaPower", "ANY", function(object, ...) 0)
setMethod("getExponentialTerm", "ANY", function(object, ...) c(0, 0))
setMethod("getPolynomialTerm", "ANY", function(object, ...) 0)
          

fixefDistributions <- c("flat", "normal", "t")
covDistributions   <- c("flat", "wishart", "invwishart",
                        "gamma", "invgamma", "custom")
residDistributions <- c("flat", "gamma", "invgamma", "point")

lmmDistributions <- list(
  flat = function() NULL,
  normal = function(sd = c(10, 2.5), cov, common.scale = TRUE) {
    matchedCall <- match.call()
    if (!is.null(matchedCall$sd)) sd <- eval(matchedCall$sd)
    if (!is.null(matchedCall$cov)) cov <- eval(matchedCall$cov)
    if (!is.null(matchedCall$sd) && !is.null(matchedCall$cov))
      warning("both sd and cov supplied to normal - only cov will be used")
    common.scale <- blme:::deparseCommonScale(common.scale)

    if (missing(cov) && !is.null(sd)) {
      sd <- sd^2
      if (length(sd) == 1) {
        cov <- diag(sd, p)
      } else if (length(sd) == 2) {
        cov <- diag(sd[c(1, rep(2, p - 1))], p)
      } else {
        sd <- rep(sd, p %/% length(sd) + 1)[1:p]
        cov <- diag(sd, p)
      }
    }
    if (missing(cov) || is.null(cov)) {
      stop("normal prior requires either sd or cov to be specified")
    }
    if (length(cov) == p) {
      cov <- diag(cov, p)
    } else if (length(cov) != p^2) {
      stop("normal prior covariance of improper length")
    }

    if (any(cov != t(cov))) stop("normal covariance not symmetric")
    
    logDet <- determinant(cov, TRUE)
    if (logDet$sign < 0 || is.infinite(logDet$modulus))
      stop("normal prior covariance negative semi-definite")
    
    new("bmerNormalDist", commonScale = common.scale, R.cov.inv = solve(chol(cov)))
  },
  t = function(df = 3, scale = c(10^2, 2.5^2), common.scale = TRUE) {
    matchedCall <- match.call()
    if (!is.null(matchedCall$df)) df <- eval(matchedCall$df)
    if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)
    common.scale <- blme:::deparseCommonScale(common.scale)

    if (df <= 0) stop("t prior requires positive degrees of freedom")

    if (length(scale) == 1) {
      scale <- diag(scale, p)
    } else if (length(scale) == 2) {
      scale <- diag(scale[c(1, rep(2, p - 1))], p)
    } else if (length(scale) == p) {
      scale <- diag(scale, p)
    } else if (length(scale) != p^2) {
      stop("t prior scale of improper length")
    }

    if (any(scale != base::t(scale))) stop("t scale not symmetric")
    
    logDet <- determinant(scale, TRUE)
    if (logDet$sign < 0 || is.infinite(logDet$modulus))
      stop("t prior scale negative semi-definite")
    
    new("bmerTDist", commonScale = common.scale, df = df, R.scale.inv = solve(chol(scale)))
  },
  gamma = function(shape = 2.5, rate = 0, common.scale = TRUE, posterior.scale = "sd") {
    matchedCall <- match.call()
    if (!is.null(matchedCall$shape)) shape <- eval(matchedCall$shape)
    if (!is.null(matchedCall$rate)) rate <- eval(matchedCall$rate)
    common.scale <- blme:::deparseCommonScale(common.scale)
    
    if (level.dim > 1) {
      warning("gamma prior applied to multivariate grouping level will be ignored")
      return(NULL)
    }

    if (shape < 0) stop("gamma requires a shape >= 0")
    if (rate  < 0) stop("gamma requires a rate >= 0")

    new("bmerGammaDist", commonScale = common.scale, shape = shape, rate = rate, posteriorScale = posterior.scale)
  },
  invgamma = function(shape = 0.001, scale = shape + 0.05, common.scale = TRUE, posterior.scale = "var") {
    matchedCall <- match.call()
    if (!is.null(matchedCall$shape)) shape <- eval(matchedCall$shape)
    if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)
    common.scale <- blme:::deparseCommonScale(common.scale)
    
    if (level.dim > 1) {
      warning("inverse gamma prior applied to multivariate grouping level will be ignored")
      return(NULL)
    }

    if (shape < 0) stop("invgamma requires a shape >= 0")
    if (scale < 0) stop("invgamma requires a scale >= 0")
    
    new("bmerInvGammaDist", commonScale = common.scale, shape = shape, scale = scale, posteriorScale = posterior.scale)
  },
  wishart = function(df = level.dim + 2.5, scale = Inf, common.scale = TRUE, posterior.scale = "cov") {
    matchedCall <- match.call()
    if (!is.null(matchedCall$df)) df <- eval(matchedCall$df)
    if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)
    common.scale <- blme:::deparseCommonScale(common.scale)
    
    if (df <= level.dim - 1)
      stop("wishart dists for degrees of freedom less than or equal to (level.dim - 1) are singular or non-existent")
    
    log.det.scale <- NULL
    if (length(scale) == 1) {
      if (is.infinite(scale)) {
        R.scale.inv <- diag(0, level.dim)
        log.det.scale <- Inf
      } else {
        if (scale[1] < 0) stop("wishart prior scale negative definite")
        R.scale.inv <- diag(1 / sqrt(scale[1]), level.dim)
      }
    } else if (length(scale) == level.dim) {
      if (any(scale < 0)) stop("wishart prior scale negative definite")
      R.scale.inv <- diag(1 / sqrt(scale), level.dim)
    } else if (length(scale) != level.dim^2) {
      stop("wishart prior scale of improper length")
    } else {
      if (all(is.infinite(scale))) {
        R.scale.inv <- diag(0, level.dim)
        log.det.scale <- Inf
      }
      R.scale.inv <- solve(chol(scale))
    }
    if (is.null(log.det.scale)) {
      if (any(diag(R.scale.inv) < 0)) stop("wishart prior scale negative definite")
      
      if (any(is.infinite(diag(R.scale.inv))))
        log.det.scale <- Inf
      else
        log.det.scale <- -2.0 * sum(log(diag(R.scale.inv)))
    }

    new("bmerWishartDist", commonScale = common.scale, df = df, R.scale.inv = R.scale.inv,
        log.det.scale = log.det.scale,
        posteriorScale = posterior.scale)
  },
  invwishart = function(df = level.dim - 0.998, scale = diag(df + 0.1, level.dim),
                        common.scale = TRUE, posterior.scale = "cov") {
    matchedCall <- match.call()
    if (!is.null(matchedCall$df)) df <- eval(matchedCall$df)
    if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)
    common.scale <- blme:::deparseCommonScale(common.scale)

    if (df <= level.dim - 1)
      stop("inverse wishart dists for degrees of freedom less than or equal to (level.dim - 1) are singular or non-existent")

    log.det.scale <- NULL
    if (length(scale) == 1) {
      if (scale == 0) {
        R.scale <- diag(0, level.dim)
        log.det.scale <- -Inf
      } else {
        if (scale[1] < 0) stop("inverse wishart prior scale negative definite")
        R.scale <- diag(sqrt(scale[1]), level.dim)
      }
    } else if (length(scale) == level.dim) {
      if (any(scale < 0)) stop("inverse wishart prior scale negative definite")
      R.scale <- diag(sqrt(scale), level.dim)
    } else if (length(scale) != level.dim^2) {
      stop("inverse wishart prior scale of improper length")
    } else {
      if (all(scale == 0)) {
        R.scale <- diag(0, level.dim)
        log.det.scale <- -Inf
      }
      R.scale <- chol(scale)
    }
    if (is.null(log.det.scale)) {
      if (any(diag(R.scale) < 0)) stop("inverse wishart prior scale negative definite")
      
      if (any(diag(R.scale) == 0))
        log.det.scale <- -Inf
      else
        log.det.scale <- 2.0 * sum(log(diag(R.scale)))
    }

    new("bmerInvWishartDist", commonScale = common.scale, df = df, R.scale = R.scale,
        log.det.scale = log.det.scale,
        posteriorScale = posterior.scale)
  },
  point = function(value = 1.0, posterior.scale = "sd") {
    matchedCall <- match.call()
    if (!is.null(matchedCall$value)) value <- eval(matchedCall$value)

    if (!(posterior.scale %in% c("sd", "var")))
      stop("point prior scale '", posterior.scale, "' unrecognized")

    if (posterior.scale == "var") value <- sqrt(value)

    if (value <= 0) stop("residual variance must be positive")
    
    new("bmerPointDist", commonScale = FALSE, value = value)
  },
  custom = function(fn, chol = FALSE, common.scale = TRUE, scale = "none") {
    matchedCall <- match.call()
    
    if (!is.null(matchedCall$chol)) chol <- eval(matchedCall$chol)
    if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)
    common.scale <- blme:::deparseCommonScale(common.scale)

    new("bmerCustomDist", fnName = matchedCall$fn, fn = fn,
        chol = chol, scale = scale, commonScale = common.scale)
  }
)

## closure out the common scale param
glmmDistributions <- list(
  flat = lmmDistributions$flat,
  normal = function(sd = c(10, 2.5), cov) {
    normal <- blme:::lmmDistributions$normal
    environment(normal) <- environment()
    
    matchedCall <- match.call()
    if (!is.null(matchedCall$sd) && !is.null(matchedCall$cov))
      warning("both sd and cov supplied to normal - only cov will be used")
    if (!is.null(matchedCall$cov)) {
      cov <- eval(matchedCall$cov)
      return(normal(cov = cov, common.scale = FALSE))
    }
    if (!is.null(matchedCall$sd)) sd <- eval(matchedCall$sd)

    normal(sd = sd, common.scale = FALSE)
  },
  t = function(df = 3, scale = c(10^2, 2.5^2)) {
    t <- blme:::lmmDistributions$t
    environment(t) <- environment()
        
    matchedCall <- match.call()
    if (!is.null(matchedCall$df)) df <- eval(matchedCall$df)
    if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)
    
    t(df = df, scale = scale, common.scale = FALSE)
  },
  gamma = function(shape = 2.5, rate = 0, posterior.scale = "sd") {
    gamma <- blme:::lmmDistributions$gamma
    environment(gamma) <- environment()
    gamma(shape, rate, TRUE, posterior.scale) 
  },
  invgamma = function(shape = 0.5, scale = 10^2, posterior.scale = "sd") {
    invgamma <- blme:::lmmDistributions$invgamma
    environment(invgamma) <- environment()
    invgamma(shape, scale, TRUE, posterior.scale)
  },
  wishart = function(df = level.dim + 2.5, scale = Inf, common.scale = TRUE, posterior.scale = "cov") {
    wishart <- blme:::lmmDistributions$wishart
    environment(wishart) <- environment()
    wishart(df, scale, TRUE, posterior.scale)
  },
  invwishart = function(df = level.dim - 0.5, scale = diag(10^2 / (df + level.dim + 1), level.dim),
                        common.scale = TRUE, posterior.scale = "cov") {
    invwishart <- blme:::lmmDistributions$invwishart
    environment(invwishart) <- environment()
    invwishart(df, scale, TRUE, posterior.scale)
  },
  custom = function(fn, chol = FALSE, scale = "none") {
    custom <- blme:::lmmDistributions$custom
    environment(custom) <- environment()
    custom(fn, chol, scale, TRUE)
  }
)

residualVarianceGammaPrior <- function(shape = 0, rate = 0, posterior.scale = "var") {
  matchedCall <- match.call()
  if (!is.null(matchedCall$shape)) shape <- eval(matchedCall$shape)
  if (!is.null(matchedCall$rate)) rate <- eval(matchedCall$rate)

  if (shape < 0) stop("gamma requires a shape >= 0")
  if (rate  < 0) stop("gamma requires a rate >= 0")
  
  new("bmerGammaDist", commonScale = FALSE, shape = shape, rate = rate, posteriorScale = posterior.scale)
}

residualVarianceInvGammaPrior <- function(shape = 0, scale = 0, posterior.scale = "var") {
  matchedCall <- match.call()
  if (!is.null(matchedCall$shape)) shape <- eval(matchedCall$shape)
  if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)

  if (shape < 0) stop("invgamma requires a shape >= 0")
  if (scale < 0) stop("invgamma requires a scale >= 0")
  
  new("bmerInvGammaDist", commonScale = FALSE, shape = shape, scale = scale, posteriorScale = posterior.scale)
}


## rather annoying problem of legacy interface allowing character strings of "true" or
## what not
deparseCommonScale <- function(common.scale) {
  if (is.null(common.scale)) return(TRUE)
  if (is.character(common.scale)) {
    if (common.scale == "TRUE" || common.scale == "true") return(TRUE)
    if (common.scale == "FALSE" || common.scale == "false") return(FALSE)
    return(eval(parse(text = common.scale)[[1]]))
  }

  common.scale
}
