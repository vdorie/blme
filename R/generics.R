setGeneric("sigma", function(object, ...) standardGeneric("sigma"))

setMethod("sigma", signature(object = "bmer"),
          function (object, ...) {
            dd <- object@dims
            if(!dd[["useSc"]]) return(1)
            object@deviance[[if(dd[["REML"]]) "sigmaREML" else "sigmaML"]]
          })


setMethod("summary", signature(object = "bmer"),
	  function(object, ...)
{
  REML <- object@dims[["REML"]]
  fcoef <- fixef(object)
  vcov <- vcov(object)
  corF <- vcov@factors$correlation
  dims <- object@dims
  coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
  llik <- logLik(object, REML)
  dev <- object@deviance
  mType <- if((non <- as.logical(length(object@V)))) "NMM" else "LMM"
  if (gen <- as.logical(length(object@muEta)))
    mType <- paste("G", mType, sep = '')
  mName <- switch(mType, LMM = "Linear", NMM = "Nonlinear",
                  GLMM = "Generalized linear",
                  GNMM = "Generalized nonlinear")
  method <- {
    if (mType == "LMM")
      if(REML) "REML" else "maximum likelihood"
    else
      paste("the", if(dims[["nAGQ"]] == 1) "Laplace" else
            "adaptive Gaussian Hermite",
            "approximation")
  }
  AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
                         logLik = as.vector(llik),
                         deviance = dev[["ML"]],
                         REMLdev = dev[["REML"]],
                         penalty = getPriorPenalty(object),
                         row.names = "")
  if (is.na(AICframe$REMLdev)) AICframe$REMLdev <- NULL
  varcor <- VarCorr(object)
  REmat <- formatVC(varcor)
  if (is.na(attr(varcor, "sc")))
    REmat <- REmat[-nrow(REmat), , drop = FALSE]
  
  if (nrow(coefs) > 0) {
    if (!dims[["useSc"]]) {
      coefs <- coefs[, 1:2, drop = FALSE]
      stat <- coefs[,1]/coefs[,2]
      pval <- 2*pnorm(abs(stat), lower.tail = FALSE)
      coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
    } else {
      stat <- coefs[,1]/coefs[,2]
      ##pval <- 2*pt(abs(stat), coefs[,3], lower.tail = FALSE)
      coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
    }
  } ## else : append columns to 0-row matrix ...
  new("summary.bmer",
      object,
      methTitle = paste(mName, "mixed model fit by", method),
      logLik = llik,
      ngrps = sapply(object@flist, function(x) length(levels(x))),
      sigma = sigma(object),
      coefs = coefs,
      vcov = vcov,
      REmat = REmat,
      AICtab= AICframe
      )
})## summary()


## This is modeled a bit after  print.summary.lm :
printBmer <- function(x, digits = max(3, getOption("digits") - 3),
                     correlation = TRUE, symbolic.cor = FALSE,
                     signif.stars = getOption("show.signif.stars"), ...)
{
  if (!inherits(x, "summary.bmer")) so <- summary(x)
  else so <- x;
  
  REML <- so@dims[["REML"]]
  llik <- so@logLik
  dev <- so@deviance
  dims <- x@dims
  
  cat(so@methTitle, "\n")
  if (!is.null(x@call$formula))
    cat("Formula:", deparse(x@call$formula),"\n")
  if (!is.null(x@call$data))
    cat("   Data:", deparse(x@call$data), "\n")
  if (!is.null(x@call$subset))
    cat(" Subset:", deparse(x@call$subset),"\n")

  covariancePriorOutput <- covariancePriorToString(x);
  if (length(covariancePriorOutput) > 0) {
    cat("Cov prior  : ", covariancePriorOutput[1], "\n", sep="");
    if (length(covariancePriorOutput) > 1) {
      for (i in 2:length(covariancePriorOutput))
        cat("           : ", covariancePriorOutput[i], "\n", sep="");
    }
  }
  unmodeledCoefficientPriorOutput <- unmodeledCoefficientPriorToString(x);
  if (length(unmodeledCoefficientPriorOutput) > 0)
    cat("Fixef prior: ", unmodeledCoefficientPriorOutput, "\n", sep="");
  
  commonScalePriorOutput <- commonScalePriorToString(x);
  if (length(commonScalePriorOutput) > 0)
    cat("Var prior : ", commonScalePriorOutput, "\n", sep="");
  
  print(so@AICtab, digits = digits)
  
  cat("Random effects:\n")
  print(so@REmat, quote = FALSE, digits = digits, ...)
  
  ngrps <- so@ngrps
  cat(sprintf("Number of obs: %d, groups: ", dims[["n"]]))
  cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
  cat("\n")
  if (is.na(so@sigma))
    cat("\nEstimated scale (compare to 1):",
        sqrt(exp(so@deviance[["lr2"]])/so@dims[["n"]]), "\n")
  if (nrow(so@coefs) > 0) {
    cat("\nFixed effects:\n")
    printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
                 digits = digits, signif.stars = signif.stars)
    if(correlation) {
      corF <- so@vcov@factors$correlation
      if (!is.null(corF)) {
        p <- ncol(corF)
        if (p > 1) {
          rn <- rownames(so@coefs)
          rns <- abbreviate(rn, minlength=11)
          cat("\nCorrelation of Fixed Effects:\n")
          if (is.logical(symbolic.cor) && symbolic.cor) {
            corf <- as(corF, "matrix")
            dimnames(corf) <- list(rns,
                                   abbreviate(rn, minlength=1, strict=TRUE))
            print(symnum(corf))
          }
          else {
            corf <- matrix(format(round(corF@x, 3), nsmall = 3),
                           ncol = p,
                           dimnames = list(rns, abbreviate(rn, minlength=6)))
            corf[!lower.tri(corf)] <- ""
            print(corf[-1, -p, drop=FALSE], quote = FALSE)
          }
        }
      }
    }
  }
  invisible(x)
}

setMethod("print", "bmer", printBmer)
setMethod("show", "bmer", function(object) printBmer(object))

formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
{
  sc <- unname(attr(varc, "sc"))
  recorr <- lapply(varc, attr, "correlation")
  reStdDev <- c(lapply(varc, attr, "stddev"), list(Residual = sc))
  reLens <- unlist(c(lapply(reStdDev, length)))
  nr <- sum(reLens)
  reMat <- array('', c(nr, 4),
                 list(rep.int('', nr),
                      c("Groups", "Name", "Variance", "Std.Dev.")))
  reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
  reMat[,2] <- c(unlist(lapply(varc, colnames)), "")
  reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
  reMat[,4] <- format(unlist(reStdDev), digits = digits)
  if (any(reLens > 1)) {
    maxlen <- max(reLens)
    corr <-
      do.call("rBind",
              lapply(recorr,
                     function(x, maxlen) {
                       x <- as(x, "matrix")
                       cc <- format(round(x, 3), nsmall = 3)
                       cc[!lower.tri(cc)] <- ""
                       nr <- dim(cc)[1]
                       if (nr >= maxlen) return(cc)
                       cbind(cc, matrix("", nr, maxlen-nr))
                     }, maxlen))
    colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
    cbind(reMat, rBind(corr, rep.int("", ncol(corr))))
  } else reMat
}

printBmerPrior <- function(x, digits = max(3, getOption("digits") - 3))
{
  hyperparameters <- x@hyperparameters;
  if (x@type == getEnumOrder(typeEnumeration, DIRECT_TYPE_NAME)) {
    if (x@families[1] == getEnumOrder(familyEnumeration, NORMAL_FAMILY_NAME)) {
      hyperparameters <- hyperparameters[-1];
    }
  }
  hyperparameters <- round(hyperparameters, digits);
  familyString <- buildStringForFamily(x@families, x@scales, hyperparameters, FALSE);
  cat(familyString$string, "\n", sep="");
}

setMethod("print", "bmerPrior", printBmerPrior);
setMethod("show", "bmerPrior", function(object) printBmerPrior(object))
