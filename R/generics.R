## This is modeled a bit after	print.summary.lm :
## Prints *both*  'mer' and 'merenv' - as it uses summary(x) mainly
printBmerenv <- function(x, digits = max(3, getOption("digits") - 3),
			correlation = NULL, symbolic.cor = FALSE,
			signif.stars = getOption("show.signif.stars"),
			ranef.comp = c("Variance", "Std.Dev."), ...)
{
  printPriors(x$priors, x$ngrps, digits)
  cat("Prior dev  : ", round(x$devcomp$cmp[["penalty"]], digits), "\n\n", sep = "")
  
  printMethod <- getS3method("print", "summary.merMod")
  result <- printMethod(x, digits, correlation, symbolic.cor, signif.stars, ranef.comp, ...)
  invisible(result)
}

print.bmerMod <- function(x, digits = max(3, getOption("digits") - 3),
                          correlation = NULL, symbolic.cor = FALSE,
                          signif.stars = getOption("show.signif.stars"),
                          ranef.comp = "Std.Dev.", ...)
{
  printPriors(x@priors, x@cnms, digits)
  cat("Prior dev  : ", round(x@devcomp$cmp[["penalty"]], digits), "\n\n", sep = "")
  
  printMethod <- getS3method("print", "merMod")
  result <- printMethod(x, digits, correlation, symbolic.cor, signif.stars, ranef.comp, ...)
  invisible(result)
}

setMethod("show",  "bmerMod", function(object) print.bmerMod(object))

print.summary.bmerMod <- printBmerenv

summary.bmerMod <- function(object, ...)
{
  result <- NextMethod(.Generic, object = object, ...)
  result$priors <- object@priors
  
  structure(result,
            class = c("summary.bmerMod", "summary.merMod"))
}

summary.summary.bmerMod <- function(object, varcov = TRUE, ...) {
  getS3method("summary", "summary.merMod")(object, varcov, ...)
}


##' @importFrom stats vcov
##' @S3method vcov summary.merMod
vcov.summary.bmerMod <- function(object, correlation = TRUE, ...) {
  getS3method("vcov", "summary.merMod")(object, correlation, ...)
}
