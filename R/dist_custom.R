setClass("bmerCustomDist",
         representation(fnName = "name",
                        fn     = "function",
                        chol   = "logical",
                        scale  = "character"),
         contains = "bmerDist",
         validity = function(object) object@scale == "log" || object@scale == "dev" || object@scale == "none")

toString.bmerCustomDist <- function(x, digits = getOption("digits"), ...) {
  paste0("custom(fn = ", x@fnName, ", ",
         "chol = ", x@chol, ", ",
         "scale = ", x@scale, ", ",
         "common.scale = ", x@commonScale, ")")
}

setMethod("getExponentialTerm", "bmerCustomDist",
  function(object, Lambda.t) {
    result <- object@fn(if (object@chol) Lambda.t else crossprod(Lambda.t))
    
    if (object@scale == "log") {
      result <- -2 * result
    } else if (object@scale == "none") {
      result <- -2 * log(result)
    }
    
    c(0, result)
  })
