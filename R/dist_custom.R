setClass("bmerCustomDist",
         slots = list(fnName = "name",
                      fn     = "function",
                      chol   = "logical",
                      scale  = "character"),
         validity = function(object) object@scale == "log" || object@scale == "dev" || object@scale == "none",
         contains = "bmerDist");

toString.bmerCustomDist <- function(x, digits = getOption("digits"), ...) {
  paste("custom(fn = ", x@fnName,
        ", chol = ", x@chol,
        ", scale = ", x@scale,
        ", common.scale = ", x@commonScale,
        ")", sep = "");
}

setMethod("getExponentialTerm", "bmerCustomDist",
  function(object, Lambda.t) {
    result <- object@fn(if (object@chol) Lambda.t else crossprod(Lambda.t))
    if (object@scale == "log") {
      result <- -2 * result
    } else if (object@scale == "none") {
      result <- -2 * log(result);
    }
    c(0, result);
  });
