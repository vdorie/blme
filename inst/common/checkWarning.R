checkWarning <- function (expr, msg = "", silent = getOption("RUnit")$silent) 
{
  tryWarn <- function (expr, silent = FALSE) {
    tryCatch(expr, warning = function(e) {
      call <- conditionCall(e)
      if (!is.null(call)) {
        if (identical(call[[1L]], quote(doTryCatch))) 
          call <- sys.call(-4L)
        dcall <- deparse(call)[1L]
        prefix <- paste("Warning in", dcall, ": ")
        LONG <- 75L
        msg <- conditionMessage(e)
        sm <- strsplit(msg, "\n")[[1L]]
        w <- 14L + nchar(dcall, type = "w") + nchar(sm[1L], 
                                  type = "w")
        if (is.na(w)) 
          w <- 14L + nchar(dcall, type = "b") + nchar(sm[1L], 
                                    type = "b")
        if (w > LONG) 
          prefix <- paste(prefix, "\n  ", sep = "")
      }
      else prefix <- "Warning : "
      msg <- paste(prefix, conditionMessage(e), "\n", sep = "")
      .Internal(seterrmessage(msg[1L]))
      if (!silent && identical(getOption("show.error.messages"), 
                               TRUE)) {
        cat(msg, file = stderr())
        .Internal(printDeferredWarnings())
      }
      invisible(structure(msg, class = "try-error"))
    })
  }

    if (missing(expr)) {
        stop("'expr' is missing")
    }
    if (is.null(silent)) {
        silent <- FALSE
        warning("'silent' has to be of type 'logical'. Was NULL. Set to FALSE.")
    }
    if (RUnit:::.existsTestLogger()) {
        .testLogger$incrementCheckNum()
    }
    if (!inherits(tryWarn(eval(expr, envir = parent.frame()), silent = silent), 
        "try-error")) {
        if (RUnit:::.existsTestLogger()) {
            .testLogger$setFailure()
        }
        stop("Warning not generated as expected\n", msg)
    }
    else {
        return(TRUE)
    }
}
