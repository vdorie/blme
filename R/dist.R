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
## note, relative to -1 / 2
if (!isGeneric("getExponentialTerm"))
  setGeneric("getExponentialTerm", function(object, ...) standardGeneric("getExponentialTerm"))
if (!isGeneric("getPolynomialTerm"))
  setGeneric("getPolynomialTerm", function(object, ...) standardGeneric("getPolynomialTerm"))

setMethod("getDFAdjustment", "ANY", function(object, ...) 0)
setMethod("getConstantTerm", "ANY", function(object, ...) 0)
setMethod("getExponentialSigmaPower", "ANY", function(object, ...) 0)
setMethod("getExponentialTerm", "ANY", function(object, ...) c(0, 0))
setMethod("getPolynomialTerm", "ANY", function(object, ...) 0)
          

fixefDistributions <- c("flat", "normal", "t", "horseshoe")
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
      if (!is.null(names(sd)) && any(names(sd) %not_in% .fixefNames & names(sd) != ""))
        stop("unrecognized fixed effects for normal prior: '",
             paste0(names(sd)[names(sd) %not_in% .fixefNames & names(sd) != ""], collapse = "', '"), "'")
      
      sd <- sd^2
      if (length(sd) == 1L) {
        cov <- if (!is.null(names(sd))) {
          diag(c(rep_len(Inf, which(.fixefNames %in% names(sd)) - 1L),
                 sd,
                 rep_len(Inf, p - which(.fixefNames %in% names(sd)))), p)
        } else { diag(sd, p) }
      } else if (length(sd) == 2L && p > 2L) {
        if (!is.null(names(sd)) && any(names(sd) != "")) {
          if (!any(names(sd) == ""))
            stop("for 2-parameter default normal prior specification with names, the default name must be left blank")
          sd <- c(rep_len(sd[names(sd) == ""], which(.fixefNames %in% names(sd)) - 1L),
                  sd[names(sd) != ""],
                  rep_len(sd[names(sd) == ""], p - which(.fixefNames %in% names(sd))))
        } else {
          sd <- c(sd[1L], rep_len(sd[2L], p - 1L))
        }
        cov <- diag(sd, p)
      } else {
        if (length(sd) > p) warning("length of sd in normal prior exceeds number of fixed effects")
        
        sd.names <- names(sd)
        
        sd <- rep_len(sd, p)
        if (!is.null(sd.names)) {
          names(sd) <- rep_len(sd.names, p)
          ind <- match(.fixefNames, names(sd))
          ind[is.na(ind)] <- which(names(sd) %not_in% .fixefNames)
          sd <- sd[ind]
        }
        cov <- diag(sd, p)
      }
    }
    if (missing(cov) || is.null(cov)) {
      stop("normal prior requires either sd or cov to be specified")
    }
    if (length(cov) == p) {
      cov <- diag(cov, p)
    } else if (length(cov) != p * p) {
      stop("normal prior covariance of improper length")
    }

    if (any(cov != t(cov))) stop("normal covariance not symmetric")
    
    logDet <- determinant(cov, TRUE)
    if (logDet$sign < 0)
      stop("normal prior covariance negative semi-definite")
    if (is.infinite(logDet$modulus)) {
      if (any(cov[upper.tri(cov) | lower.tri(cov)] != 0))
        stop("normal prior covariance infinite")
      ## special case for diagonal scale matrices with infinite variances
      R.cov.inv <- diag(1 / sqrt(diag(cov)))
    } else {
      R.cov.inv <- solve(chol(cov))
    }
    
    new("bmerNormalDist", commonScale = common.scale, R.cov.inv = R.cov.inv)
  },
  t = function(df = 3, mean = 0, scale = c(10^2, 2.5^2), common.scale = TRUE) {
    matchedCall <- match.call()
    if (!is.null(matchedCall$df)) df <- eval(matchedCall$df)
    if (!is.null(matchedCall$mean)) mean <- eval(matchedCall$mean)
    if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)
    common.scale <- blme:::deparseCommonScale(common.scale)

    if (df <= 0) stop("t prior requires positive degrees of freedom")
    
    if (length(mean) == 1L) {
      mean <- rep_len(mean, p)
    } else if (length(mean) == 2L && p > 1L) {
      mean <- c(mean[1L], rep_len(mean[2L], p - 1L))
    } else if (length(mean) != p) {
      stop("t prior mean of improper length")
    }
    
    if (length(scale) == 1L) {
      scale <- diag(scale, p)
    } else if (length(scale) == 2L && p > 1L) {
      scale <- diag(c(scale[1L], rep_len(scale[2L], p - 1L)), p)
    } else if (length(scale) == p) {
      scale <- diag(scale, p)
    } else if (length(scale) != p * p) {
      stop("t prior scale of improper length")
    }
    
    if (any(scale != base::t(scale))) stop("t scale not symmetric")
    
    logDet <- determinant(scale, TRUE)
    if (logDet$sign < 0)
      stop("t prior scale negative semi-definite")
    if (is.infinite(logDet$modulus)) {
      if (any(scale[upper.tri(scale) | lower.tri(scale)] != 0))
        stop("t prior scale infinite")
      ## special case for diagonal scale matrices with infinite variances
      R.scale.inv <- diag(1 / sqrt(diag(scale)))
      d <- sum(is.finite(diag(scale)))
    } else {
      R.scale.inv <- solve(chol(scale))
      d <- nrow(R.scale.inv)
    }
    
    new("bmerTDist", commonScale = common.scale, df = df, beta.0 = mean, d = d, R.scale.inv = R.scale.inv)
  },
  horseshoe = function(mean = 0, global.shrinkage = 2.5, common.scale = TRUE) {
    matchedCall <- match.call()
    if (!is.null(matchedCall$mean)) mean <- eval(matchedCall$mean)
    if (!is.null(matchedCall$global.shrinkage)) global.shrinkage <- eval(matchedCall$global.shrinkage)
    common.scale <- blme:::deparseCommonScale(common.scale)
    
    if (length(mean) == 1L) {
      mean <- rep_len(mean, p)
    } else if (length(mean) == 2L && p > 1L) {
      mean <- c(mean[1L], rep_len(mean[2L], p - 1L))
    } else if (length(mean) != p) {
      stop("horseshoe prior mean of improper length")
    }
    
    if (global.shrinkage <= 0)
      stop("horseshoe prior global.shrinkage parameter must be positive")
    if (length(global.shrinkage) != 1L)
      stop("horseshoe prior global.shrinkage must be of length 1")
    
    new("bmerHorseshoeDist", beta.0 = mean, tau.sq = global.shrinkage^2, commonScale = common.scale)
  },
  gamma = function(shape = 2.5, rate = 0, common.scale = TRUE, posterior.scale = "sd") {
    matchedCall <- match.call()
    if (!is.null(matchedCall$shape)) shape <- eval(matchedCall$shape)
    if (!is.null(matchedCall$rate)) rate <- eval(matchedCall$rate)
    common.scale <- blme:::deparseCommonScale(common.scale)
    
    if (level.dim > 1L) {
      warning("gamma prior applied to multivariate grouping level will be ignored")
      return(NULL)
    }

    if (shape < 0) stop("gamma prior shape must be positive")
    if (rate  < 0) stop("gamma prior rate must be positive")

    new("bmerGammaDist", commonScale = common.scale, shape = shape, rate = rate, posteriorScale = posterior.scale)
  },
  invgamma = function(shape = 0.001, scale = shape + 0.05, common.scale = TRUE, posterior.scale = "var") {
    matchedCall <- match.call()
    if (!is.null(matchedCall$shape)) shape <- eval(matchedCall$shape)
    if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)
    common.scale <- blme:::deparseCommonScale(common.scale)
    
    if (level.dim > 1L) {
      warning("inverse gamma prior applied to multivariate grouping level will be ignored")
      return(NULL)
    }

    if (shape < 0) stop("invgamma prior shape must be positive")
    if (scale < 0) stop("invgamma prior scale must be positive")
    
    new("bmerInvGammaDist", commonScale = common.scale, shape = shape, scale = scale, posteriorScale = posterior.scale)
  },
  wishart = function(df = level.dim + 2.5, scale = Inf, common.scale = TRUE, posterior.scale = "cov") {
    matchedCall <- match.call()
    if (!is.null(matchedCall$df)) df <- eval(matchedCall$df)
    if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)
    common.scale <- blme:::deparseCommonScale(common.scale)
    
    if (df <= level.dim - 1L)
      stop("wishart dists for degrees of freedom less than or equal to (level.dim - 1) are singular or non-existent")
    
    log.det.scale <- NULL
    if (length(scale) == 1L) {
      if (is.infinite(scale)) {
        R.scale.inv <- diag(0, level.dim)
        log.det.scale <- Inf
      } else {
        if (scale[1L] < 0) stop("wishart prior scale negative definite")
        R.scale.inv <- diag(1 / sqrt(scale[1L]), level.dim)
      }
    } else if (length(scale) == level.dim) {
      if (any(scale < 0)) stop("wishart prior scale negative definite")
      R.scale.inv <- diag(1 / sqrt(scale), level.dim)
    } else if (length(scale) != level.dim * level.dim) {
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

    if (df <= level.dim - 1L)
      stop("inverse wishart dists for degrees of freedom less than or equal to (level.dim - 1) are singular or non-existent")

    log.det.scale <- NULL
    if (length(scale) == 1L) {
      if (scale == 0) {
        R.scale <- diag(0, level.dim)
        log.det.scale <- -Inf
      } else {
        if (scale[1L] < 0) stop("inverse wishart prior scale negative definite")
        R.scale <- diag(sqrt(scale[1L]), level.dim)
      }
    } else if (length(scale) == level.dim) {
      if (any(scale < 0)) stop("inverse wishart prior scale negative definite")
      R.scale <- diag(sqrt(scale), level.dim)
    } else if (length(scale) != level.dim * level.dim) {
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
    .prior <- blme:::lmmDistributions$normal
    environment(.prior) <- environment()
    
    matchedCall <- match.call()
    if (!is.null(matchedCall$sd) && !is.null(matchedCall$cov))
      warning("both sd and cov supplied to normal - only cov will be used")
    if (!is.null(matchedCall$cov)) {
      eval(substitute(.prior(cov = cov, common.scale = FALSE), environment()),
           parent.frame())
    } else {
      eval(substitute(.prior(sd = sd, common.scale = FALSE), environment()),
           parent.frame())
    }
  },
  t = function(df = 3, mean = 0, scale = c(10^2, 2.5^2)) {
    .prior <- blme:::lmmDistributions$t
    environment(.prior) <- environment()
    
    eval(substitute(.prior(df, mean, scale, FALSE), environment()),
         parent.frame()) 
  },
  horseshoe = function(mean = 0, global.shrinkage = 2.5) {
    .prior <- blme:::lmmDistributions$horseshoe
    environment(.prior) <- environment()
    
    eval(substitute(.prior(mean, global.shrinkage, FALSE), environment()),
         parent.frame()) 
  },
  gamma = function(shape = 2.5, rate = 0, posterior.scale = "sd") {
    .prior <- blme:::lmmDistributions$gamma
    environment(.prior) <- environment()
    
    eval(substitute(.prior(shape, rate, TRUE, posterior.scale), environment()),
         parent.frame()) 
  },
  invgamma = function(shape = 0.5, scale = 10^2, posterior.scale = "sd") {
    .prior <- blme:::lmmDistributions$invgamma
    environment(.prior) <- environment()
    eval(substitute(.prior(shape, scale, TRUE, posterior.scale), environment()),
         parent.frame())
  },
  wishart = function(df = level.dim + 2.5, scale = Inf, posterior.scale = "cov") {
    .prior <- blme:::lmmDistributions$wishart
    environment(.prior) <- environment()
    
    eval(substitute(.prior(df, scale, TRUE, posterior.scale), environment()),
         parent.frame())
  },
  invwishart = function(df = level.dim - 0.5, scale = diag(10^2 / (df + level.dim + 1), level.dim),
                        posterior.scale = "cov") {
    .prior <- blme:::lmmDistributions$invwishart
    environment(.prior) <- environment()
    eval(substitute(.prior(df, scale, common.scale, posterior.scale), environment()),
         parent.frame())
  },
  custom = function(fn, chol = FALSE, scale = "none") {
    .prior <- blme:::lmmDistributions$custom
    environment(.prior) <- environment()
    eval(substitute(.prior(fn, chol, scale), environment()),
         parent.frame())
  }
)

residualVarianceGammaPrior <- function(shape = 0, rate = 0, posterior.scale = "var") {
  matchedCall <- match.call()
  if (!is.null(matchedCall$shape)) shape <- eval(matchedCall$shape)
  if (!is.null(matchedCall$rate)) rate <- eval(matchedCall$rate)

  if (shape < 0) stop("gamma prior shape must be positive")
  if (rate  < 0) stop("gamma prior rate must be positive")
  
  new("bmerGammaDist", commonScale = FALSE, shape = shape, rate = rate, posteriorScale = posterior.scale)
}

residualVarianceInvGammaPrior <- function(shape = 0, scale = 0, posterior.scale = "var") {
  matchedCall <- match.call()
  if (!is.null(matchedCall$shape)) shape <- eval(matchedCall$shape)
  if (!is.null(matchedCall$scale)) scale <- eval(matchedCall$scale)

  if (shape < 0) stop("invgamma prior shape must be positive")
  if (scale < 0) stop("invgamma prior scale must be positive")
  
  new("bmerInvGammaDist", commonScale = FALSE, shape = shape, scale = scale, posteriorScale = posterior.scale)
}


## rather annoying problem of legacy interface allowing character strings of "true" or
## what not
deparseCommonScale <- function(common.scale) {
  if (is.null(common.scale)) return(TRUE)
  if (is.character(common.scale)) {
    if (common.scale == "TRUE" || common.scale == "true") return(TRUE)
    if (common.scale == "FALSE" || common.scale == "false") return(FALSE)
    return(eval(parse(text = common.scale)[[1L]]))
  }

  common.scale
}
