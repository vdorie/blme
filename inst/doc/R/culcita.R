culcitaPath <- file.path(dataDir, "culcita.RData")
if (!file.exists(culcitaPath))
  download.file("http://glmm.wdfiles.com/local--files/examples/culcita.RData", culcitaPath)

load(culcitaPath)
rm(culcitaPath)

culcita <- culcita_dat; rm(culcita_dat)

## applies to a call to model.matrix
standardizeDesign <- function(X)
{
  standardize.binary <- function(x) x - mean(x, na.rm = TRUE)
  standardize.cont <- function(x) 0.5 * (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  standardize <- function(x) {
    y <- unique(x)
    y <- y[!is.na(y)]
    if (length(y) == 2) return(standardize.binary(x))
    standardize.cont(x)
  }
  
  X <- X[,colnames(X) != "(Intercept)"]
  apply(X, 2, standardize)
}

culcita.z <- with(culcita,
  data.frame(block, predation, standardizeDesign(model.matrix(~ttt, culcita))))

culcitaSep <- culcita[-c(19, 20),]
culcitaSep.z <- with(culcitaSep,
  data.frame(block, predation, standardizeDesign(model.matrix(~ttt, culcitaSep))))

culcitaFitPath <- file.path(dataDir, "culcitaFit.RData")
if (!file.exists(culcitaFitPath)) {
  m1 <- glmer(predation ~ tttcrabs + tttshrimp + tttboth + (1 | block),
              culcita.z, family = binomial, nAGQ = 10)
  
  m2 <- glmer(predation ~ tttcrabs + tttshrimp + tttboth + (1 | block),
              culcitaSep.z, family = binomial, nAGQ = 10)
  
  m3 <- bglmer(predation ~ tttcrabs + tttshrimp + tttboth + (1 | block),
               culcitaSep.z, family = binomial, nAGQ = 10,
               cov.prior = NULL, fixef.prior = t(1, 0.75))

  save(m1, m2, m3, file = culcitaFitPath)
} else load(culcitaFitPath)

rm(culcitaFitPath, standardizeDesign)

invlogit <- function(x) { e.x <- exp(x); e.x / (1 + e.x) }


imgPath <- file.path(imgDir, "culcita.pdf")
if (!file.exists(imgPath)) {
  m1DevFun <- glmer(predation ~ tttcrabs + tttshrimp + tttboth + (1 | block),
                    culcita.z, family = binomial, devFunOnly = TRUE, nAGQ = 10)
  body(m1DevFun) <- expression({
    resp$setOffset(baseOffset)
    resp$updateMu(lp0)
    pp$setTheta(as.double(theta))
    spars <- as.numeric(pars[-dpars])
    offset <- if (length(spars) == 0) 
      baseOffset
    else baseOffset + pp$X %*% spars
    resp$setOffset(offset)
    p <- pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, 
                     verbose)
    resp$updateWts()
    p
  })
  environment(m1DevFun)$fixef <- m1@beta


  m2DevFun <- glmer(predation ~ tttcrabs + tttshrimp + tttboth + (1 | block),
                    culcitaSep.z, family = binomial, devFunOnly = TRUE, nAGQ = 10)
  body(m2DevFun) <- body(m1DevFun)
  environment(m2DevFun)$fixef <- m2@beta

  m3DevFun <- bglmer(predation ~ tttcrabs + tttshrimp + tttboth + (1 | block),
                     culcitaSep.z, family = binomial, devFunOnly = TRUE, nAGQ = 10,
                     cov.prior = NULL, fixef.prior = t(1, 0.75))
  body(m3DevFun) <- expression({
    resp$setOffset(baseOffset)
    resp$updateMu(lp0)
    pp$setTheta(as.double(theta))
    spars <- as.numeric(pars[-dpars])
    offset <- if (length(spars) == 0) 
        baseOffset
    else baseOffset + pp$X %*% spars
    resp$setOffset(offset)
    p <- pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, 
        verbose)
    resp$updateWts()
    Lambda.ts <- getCovBlocks(pp$Lambdat, ranefStructure)
    exponentialTerms <- calculatePriorExponentialTerms(priors, 
        spars, Lambda.ts)
    polynomialTerm <- calculatePriorPolynomialTerm(priors$covPriors, 
        Lambda.ts)
    p + exponentialTerms[[1]] + polynomialTerm + blmerControl$constant
  })
  environment(m3DevFun)$fixef <- m3@beta


  expNormalize <- function(x, del) {
    if (length(x) == 1) return(x)
    x <- exp(x - median(x))
    x / sum(x * del)
  }
  curveWrapper <- function(x) {
    result <- sapply(x, function(x.i) {
      rho <- environment(devFun)
      rho$theta <- x.i
      opt <- optimizeGlmer(devFun, optimizer = "Nelder_Mead", 
                           restart_edge = FALSE, control = list(), 
                           start = list(theta = x.i, fixef = rho$fixef), nAGQ = 1, verbose = 0L, stage = 2)
      rho$fixef <- opt$par[-1]
      -0.5 * opt$fval
    })
    expNormalize(result, c(x[2] - x[1], diff(x)))
  }
  
  
  
  xValues <- seq(0.05, 15, length.out = 101)
  devFun <- m1DevFun; yValues1 <- curveWrapper(xValues)
  devFun <- m2DevFun; yValues2 <- curveWrapper(xValues) / 60
  devFun <- m3DevFun; yValues3 <- curveWrapper(xValues)

  pdf(imgPath, defaultImgWidth, defaultImgHeight)
  par(defaultPars)
  plot(NULL, type = "n", xlim = c(0, max(xValues)), ylim = range(yValues1, yValues2, yValues3),
       main = "Profiled Objective Fns", ylab = "density", xlab = expression(sigma[u]),
       yaxt = "n", bty = "n", yaxs = "i")
  lines(xValues, yValues1)
  lines(xValues, yValues2, lty = 2)
  lines(xValues, yValues3, col = "gray")

  legend("topright", c("lik-cmpl data", "lik-miss obs", "post-miss obs"),
         lty = c(1, 2, 1), col = c("black", "black", "gray"), bty = "n",
         cex = 0.7)
  dev.off()

  rm(xValues, yValues1, yValues2, yValues3,
     curveWrapper, expNormalize,
     devFun, m1DevFun, m2DevFun, m3DevFun)
}
rm(imgPath)
