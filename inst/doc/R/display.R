## from the package arm by Yu-Sung Su
## http://cran.r-project.org/web/packages/arm/index.html
## published under GPL v2

if (!isGeneric("display")) {
    setGeneric("display",
               function(object, ...)
               standardGeneric("display"))
}

fround <- function (x, digits) {
    format (round (x, digits), nsmall=digits)
}
  
pfround <- function (x, digits) {
    print (fround (x, digits), quote=FALSE)
}

as.matrix.VarCorr <- function (varc, useScale, digits){
# VarCorr function for lmer objects, altered as follows:
#   1.  specify rounding
#   2.  print statement at end is removed
#   3.  reMat is returned
#   4.  last line kept in reMat even when there's no error term
                  sc <- attr(varc, "sc")[[1]]   
                  if(is.na(sc)) sc <- 1 
#                  recorr <- lapply(varc, function(el) el@factors$correlation)
                  recorr <- lapply(varc, function(el) attr(el, "correlation"))
                  #reStdDev <- c(lapply(recorr, slot, "sd"), list(Residual = sc))
                  reStdDev <- c(lapply(varc, function(el) attr(el, "stddev")), list(Residual = sc))
                  reLens <- unlist(c(lapply(reStdDev, length)))
                  reMat <- array('', c(sum(reLens), 4),
                                 list(rep('', sum(reLens)),
                                      c("Groups", "Name", "Variance", "Std.Dev.")))
                  reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
                  reMat[,2] <- c(unlist(lapply(reStdDev, names)), "")
#                  reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
#                  reMat[,4] <- format(unlist(reStdDev), digits = digits)
                  reMat[,3] <- fround(unlist(reStdDev)^2, digits)
                  reMat[,4] <- fround(unlist(reStdDev), digits)
                  if (any(reLens > 1)) {
                      maxlen <- max(reLens)
                      corr <-
                          do.call("rbind",
                                  lapply(recorr,
                                         function(x, maxlen) {
                                             x <- as(x, "matrix")
#                                             cc <- format(round(x, 3), nsmall = 3)
                                             cc <- fround (x, digits)
                                             cc[!lower.tri(cc)] <- ""
                                             nr <- dim(cc)[1]
                                             if (nr >= maxlen) return(cc)
                                             cbind(cc, matrix("", nr, maxlen-nr))
                                         }, maxlen))
                      colnames(corr) <- c("Corr", rep("", maxlen - 1))
                      reMat <- cbind(reMat, rbind(corr, rep("", ncol(corr))))
                  }
#                  if (!useScale) reMat <- reMat[-nrow(reMat),]
          if (useScale<0) reMat[nrow(reMat),] <- c ("No residual sd", rep("",ncol(reMat)-1))
          return (reMat)
      }



setMethod("display", signature(object = "merMod"),
    function(object, digits=2, detail=FALSE)
    {
    out <- NULL
    out$call <- object@call
    print (out$call)
    #object <- summary(object)
    #summ <- summary(object)
    fcoef <- fixef(object)
    #coefs <- attr(summ, "coefs")
    #useScale <- attr (VarCorr (object), "sc")
    useScale <- getME(object, "devcomp")$dims["useSc"]
    corF <- vcov(object)@factors$correlation
    coefs <- cbind(fcoef, corF@sd)
    if (length (fcoef) > 0){
      if (!useScale) {
        coefs <- coefs[, 1:2, drop = FALSE]
        out$z.value <- coefs[, 1]/coefs[, 2]
        out$p.value <- 2 * pnorm(abs(out$z.value), lower = FALSE)
        coefs <- cbind(coefs, `z value` = out$z.value, `Pr(>|z|)` = out$p.value)
      }
      else {
        out$t.value <- coefs[, 1]/coefs[, 2]
        coefs <- cbind(coefs, `t value` = out$t.value)
      }
    dimnames(coefs)[[2]][1:2] <- c("coef.est", "coef.se")
      if(detail){
        pfround (coefs, digits)
      }
      else{
        pfround(coefs[,1:2], digits)
      }
    }
    out$coef <- coefs[,"coef.est"]
    out$se <- coefs[,"coef.se"]
    cat("\nError terms:\n")
    vc <- as.matrix.VarCorr (VarCorr (object), useScale=useScale, digits)
    print (vc[,c(1:2,4:ncol(vc))], quote=FALSE)
    out$ngrps <- lapply(object@flist, function(x) length(levels(x)))
    is_REML <- isREML(object)
    llik <- logLik(object, REML=is_REML)
    out$AIC <- AIC(llik)
    out$deviance <- deviance(refitML(object))     # Dbar
    out$n <- getME(object, "devcomp")$dims["n"]
    Dhat <- -2*(llik) # Dhat
    pD <- out$deviance - Dhat              # pD
    out$DIC <- out$deviance + pD               # DIC=Dbar+pD=Dhat+2pD
    cat("---\n")
    cat(sprintf("number of obs: %d, groups: ", out$n))
    cat(paste(paste(names(out$ngrps), out$ngrps, sep = ", "), collapse = "; "))
    cat(sprintf("\nAIC = %g, DIC = ", round(out$AIC,1)))
    cat(round(out$DIC, 1))
    cat("\ndeviance =", fround (out$deviance, 1), "\n")
    if (useScale < 0){
      out$sigma.hat <- .Call("mer_sigma", object, FALSE, PACKAGE = "lme4")
      cat("overdispersion parameter =", fround (out$sigma.hat, 1), "\n")
    }
    return(invisible(out))
    }
)

