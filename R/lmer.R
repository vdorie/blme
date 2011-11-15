# lmer, glmer and nlmer plus methods and utilities


### Utilities for parsing the mixed model formula

findbars <- function(term)
### Return the pairs of expressions that separated by vertical bars
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

nobars <- function(term)
### Return the formula omitting the pairs of expressions that are
### separated by vertical bars
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
	nb <- nobars(term[[2]])
	if (is.null(nb)) return(NULL)
	term[[2]] <- nb
	return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

subbars <- function(term)
### Substitute the '+' function for the '|' function
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
	term[[2]] <- subbars(term[[2]])
	return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
	term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

subnms <- function(term, nlist)
### Substitute any names from nlist in term with 1
{
    if (!is.language(term)) return(term)
    if (is.name(term)) {
        if (any(unlist(lapply(nlist, get("=="), term)))) return(1)
        return(term)
    }
    stopifnot(length(term) >= 2)
    for (j in 2:length(term)) term[[j]] <- subnms(term[[j]], nlist)
    term
}

slashTerms <- function(x)
### Return the list of '/'-separated terms in an expression that
### contains slashes
{
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

makeInteraction <- function(x)
### from a list of length 2 return recursive interaction terms
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

expandSlash <- function(bb)
### expand any slashes in the grouping factors returned by findbars
{
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}

### Utilities used in lmer, glmer and nlmer

createCm <- function(A, s)
### Create the nonzero pattern for the sparse matrix Cm from A.
### ncol(A) is s * ncol(Cm).  The s groups of ncol(Cm) consecutive
### columns in A are overlaid to produce Cm.
{
    stopifnot(is(A, "dgCMatrix"))
    s <- as.integer(s)[1]
    if (s == 1L) return(A)
    if ((nc <- ncol(A)) %% s)
        stop(gettextf("ncol(A) = %d is not a multiple of s = %d",
                      nc, s))
    ncC <- as.integer(nc / s)
    TA <- as(A, "TsparseMatrix")
    as(new("dgTMatrix", Dim = c(nrow(A), ncC),
           i = TA@i, j = as.integer(TA@j %% ncC), x = TA@x),
       "CsparseMatrix")
}

### FIXME: somehow the environment of the mf formula does not have
### .globalEnv in its parent list.  example(Mmmec, package = "mlmRev")
### used to have a formula of ~ offset(log(expected)) + ... and the
### offset function was not found in eval(mf, parent.frame(2))
lmerFrames <- function(mc, formula, contrasts, vnms = character(0))
### Create the model frame, X, Y, wts, offset and terms

### mc - matched call of calling function
### formula - two-sided formula
### contrasts - contrasts argument
### vnms - names of variables to be included in the model frame
{
    mf <- mc
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]

    ## The model formula for evaluation of the model frame.  It looks
    ## like a linear model formula but includes any random effects
    ## terms and any names of parameters used in a nonlinear mixed model.
    frame.form <- subbars(formula)      # substitute `+' for `|'
    if (length(vnms) > 0)               # add the variables names for nlmer
        frame.form[[3]] <-
            substitute(foo + bar,
                       list(foo = parse(text = paste(vnms, collapse = ' + '))[[1]],
                            bar = frame.form[[3]]))

    ## The model formula for the fixed-effects terms only.
    fixed.form <- nobars(formula)       # remove any terms with `|'
    if (inherits(fixed.form, "name"))   # RHS is empty - use `y ~ 1'
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))

    ## attach the correct environment
    environment(fixed.form) <- environment(frame.form) <- environment(formula)

    ## evaluate a model frame
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fe <- mf                            # save a copy of the call
    mf <- eval(mf, parent.frame(2))

    ## evaluate the terms for the fixed-effects only (used in anova)
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2)) # allow model.frame to update them

    ## response vector
    Y <- model.response(mf, "any")
    ## avoid problems with 1D arrays, but keep names
    if(length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
    }
    mt <- attr(fe, "terms")

    ## Extract X checking for a null model. This check shouldn't be
    ## needed because an empty formula is changed to ~ 1 but it can't hurt.
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)
    storage.mode(X) <- "double"      # when ncol(X) == 0, X is logical
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL

    ## Extract the weights and offset.  For S4 classes we want the
    ## `not used' condition to be numeric(0) instead of NULL
    wts <- model.weights(mf); if (is.null(wts)) wts <- numeric(0)
    off <- model.offset(mf); if (is.null(off)) off <- numeric(0)

    ## check weights and offset
    if (any(wts <= 0))
        stop(gettextf("negative weights or weights of zero are not allowed"))
    if(length(off) && length(off) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(off), NROW(Y)))

    ## remove the terms attribute from mf
    attr(mf, "terms") <- mt
    list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), mf = mf, fixef = fixef)
}

##' Is f1 nested within f2?
##'
##' Does every level of f1 occur in conjunction with exactly one level
##' of f2? The function is based on converting a triplet sparse matrix
##' to a compressed column-oriented form in which the nesting can be
##' quickly evaluated.
##'
##' @param f1 factor 1
##' @param f2 factor 2

##' @return TRUE if factor 1 is nested within factor 2

isNested <- function(f1, f2)
{
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    stopifnot(length(f1) == length(f2))
    sm <- as(new("ngTMatrix",
                 i = as.integer(f2) - 1L,
                 j = as.integer(f1) - 1L,
                 Dim = c(length(levels(f2)),
                 length(levels(f1)))),
             "CsparseMatrix")
    all(diff(sm@p) < 2)
}

##' dimsNames and devNames are in the package's namespace rather than
##' in the function lmerFactorList because the function sparseRasch
##' needs to access them.

dimsNames <- c("nt", "n", "p", "q", "s", "np", "LMM", "REML",
               "fTyp", "lTyp", "vTyp", "nest", "useSc", "nAGQ",
               "verb", "mxit", "mxfn", "cvg")
dimsDefault <- list(s = 1L,             # identity mechanistic model
                    mxit= 300L,         # maximum number of iterations
                    mxfn= 900L, # maximum number of function evaluations
                    verb= 0L,           # no verbose output
                    np= 0L,             # number of parameters in ST
                    LMM= 0L,            # not a linear mixed model
                    REML= 0L,         # glmer and nlmer don't use REML
                    fTyp= 2L,           # default family is "gaussian"
                    lTyp= 5L,           # default link is "identity"
                    vTyp= 1L, # default variance function is "constant"
                    useSc= 1L, # default is to use the scale parameter
                    nAGQ= 1L,                  # default is Laplace
                    cvg = 0L)                  # no optimization yet attempted

devNames <- c("ML", "REML", "ldL2", "ldRX2", "sigmaML",
              "sigmaREML", "pwrss", "disc", "usqr", "wrss",
              "dev", "llik", "NULLdev")


##' Create model matrices from r.e. terms.
##'
##' Create the list of model matrices from the random-effects terms in
##' the formula and the model frame.
##'
##' @param formula model formula
##' @param fr: list with '$mf': model frame; '$X': .. matrix
##' @param rmInt logical scalar - should the `(Intercept)` column
##'        be removed before creating Zt
##' @param drop logical scalar indicating if elements with numeric
##'        value 0 should be dropped from the sparse model matrices
##'
##' @return a list with components named \code{"trms"}, \code{"fl"}
##'        and \code{"dims"}
lmerFactorList <- function(formula, fr, rmInt, drop)
{
    mf <- fr$mf
    ## record dimensions and algorithm settings

    ## create factor list for the random effects
    bars <- expandSlash(findbars(formula[[3]]))
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars,
                 function(x)
                 {
                   ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                                         list(fac = x[[3]])), mf)
                   im <- as(ff, "sparseMatrix") # transpose of indicators
                   ## Could well be that we should rather check earlier .. :
                   if(!isTRUE(validObject(im, test=TRUE)))
		     stop("invalid conditioning factor in random effect: ", format(x[[3]]))
                   
                   mm <- model.matrix(eval(substitute(~ expr, # model matrix
                                                      list(expr = x[[2]]))),
                                      mf)
                   if (rmInt) {
                     if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
                     if (ncol(mm) < 2)
                       stop("lhs of a random-effects term cannot be an intercept only")
                     mm <- mm[ , -icol , drop = FALSE]
                   }
                   ans <- list(f = ff,
                               A = do.call(rBind,
                                 lapply(seq_len(ncol(mm)), function(j) im)),
                               Zt = do.call(rBind,
                                 lapply(seq_len(ncol(mm)),
                                        function(j) {im@x <- mm[,j]; im})),
                               ST = matrix(0, ncol(mm), ncol(mm),
                                 dimnames = list(colnames(mm), colnames(mm))))
                   if (drop) {
                     ## This is only used for nlmer models.
                     ## Need to do something more complicated for A
                     ## here.  Essentially you need to create a copy
                     ## of im for each column of mm, im@x <- mm[,j],
                     ## create the appropriate number of copies,
                     ## prepend matrices of zeros, then rBind and drop0.
                     ans$A@x <- rep(0, length(ans$A@x))
                     ans$Zt <- drop0(ans$Zt)
                   }
                   ans
                 })
    dd <-
      VecFromNames(dimsNames, "integer",
                     c(list(n = nrow(mf), p = ncol(fr$X), nt = length(fl),
                            q = sum(sapply(fl, function(el) nrow(el$Zt)))),
                       dimsDefault))
    ## order terms by decreasing number of levels in the factor but don't
    ## change the order if this is already true
    nlev <- sapply(fl, function(el) length(levels(el$f)))
    ## determine the number of random effects at this point
    if (any(diff(nlev)) > 0) fl <- fl[rev(order(nlev))]
    ## separate the terms from the factor list
    trms <- lapply(fl, "[", -1)
    names(trms) <- NULL
    fl <- lapply(fl, "[[", "f")
    attr(fl, "assign") <- seq_along(fl)
    ## check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        ## check that the lengths of the number of levels coincide
        fl <- fl[match(ufn, fnms)]
        attr(fl, "assign") <- match(fnms, ufn)
    }
    names(fl) <- ufn
    ## check for nesting of factors
    dd["nest"] <- all(sapply(seq_along(fl)[-1],
                             function(i) isNested(fl[[i-1]], fl[[i]])))

    list(trms = trms, fl = fl, dims = dd)
}

checkSTform <- function(ST, STnew)
### Check that the 'STnew' argument matches the form of ST.
{
  stopifnot(is.list(STnew), length(STnew) == length(ST),
            all.equal(names(ST), names(STnew)))
  lapply(seq_along(STnew), function (i)
         stopifnot(class(STnew[[i]]) == class(ST[[i]]),
                   all.equal(dim(STnew[[i]]), dim(ST[[i]]))))
  all(unlist(lapply(STnew, function(m) all(diag(m) > 0))))
}

lmerControl <- function(msVerbose = getOption("verbose"),
                        maxIter = 300L, maxFN = 900L)
### Control parameters for lmer, glmer and nlmer
{
    stopifnot(maxIter >= 0, maxFN >= 0)
    list(
         maxIter = as.integer(maxIter),
         maxFN = as.integer(maxFN),
	 msVerbose = as.integer(msVerbose))# "integer" on purpose
}

VecFromNames <- function(nms, mode = "numeric", defaults = list())
### Generate a named vector of the given mode
{
    ans <- vector(mode = mode, length = length(nms))
    names(ans) <- nms
    ans[] <- NA
    if ((nd <- length(defaults <- as.list(defaults))) > 0) {
        if (length(dnms <- names(defaults)) < nd)
            stop("defaults must be a named list")
        stopifnot(all(dnms %in% nms))
        ans[dnms] <- as(unlist(defaults), mode)
    }
    ans
}

mkZt <- function(FL, start, s = 1L)
### Create the standard versions of flist, Zt, Gp, ST, A, Cm,
### Cx, and L. Update dd.
{
    dd <- FL$dims
    fl <- FL$fl
    asgn <- attr(fl, "assign")
    trms <- FL$trms
    ST <- lapply(trms, `[[`, "ST")
    names(ST) <- names(fl)[attr(fl, "assign")];
    Ztl <- lapply(trms, `[[`, "Zt")
    Zt <- do.call(rBind, Ztl)
    Zt@Dimnames <- vector("list", 2)
    Gp <- unname(c(0L, cumsum(sapply(Ztl, nrow))))
    .Call(mer_ST_initialize, ST, Gp, Zt)
    A <- do.call(rBind, lapply(trms, `[[`, "A"))
    rm(Ztl, FL)                         # because they could be large
    nc <- sapply(ST, ncol)         # of columns in els of ST
    Cm <- createCm(A, s)
    L <- .Call(mer_create_L, Cm)
    if (s < 2) Cm <- new("dgCMatrix")
    if (!is.null(start) && checkSTform(ST, start)) ST <- start

    nvc <- sapply(nc, function (qi) (qi * (qi + 1))/2) # no. of var. comp.
### FIXME: Check number of variance components versus number of
### levels in the factor for each term. Warn or stop as appropriate

    dd["np"] <- as.integer(sum(nvc))    # number of parameters in optimization
    dev <- VecFromNames(devNames, "numeric")
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn

    list(Gp = Gp, ST = ST, A = A, Cm = Cm, L = L, Zt = Zt,
         dd = dd, dev = dev, flist = fl)
}

famNms <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
            "poisson", "quasibinomial", "quasipoisson", "quasi")
linkNms <- c("logit", "probit", "cauchit", "cloglog", "identity",
             "log", "sqrt", "1/mu^2", "inverse")
varNms <- c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")

famType <- function(family)
{
    if (!(fTyp <- match(family$family, famNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    if (!(lTyp <- match(family$link, linkNms, nomatch = 0)))
        stop(gettextf("unknown link: %s",
                      sQuote(family$link), domain = "R-lme4"))
    vNam <- switch(fTyp,
                   "mu(1-mu)",          # binomial
                   "constant",          # gaussian
                   "mu^2",              # Gamma
                   "mu^3",              # inverse.gaussian
                   "mu",                # poisson
                   "mu(1-mu)",          # quasibinomial
                   "mu",                # quasipoisson
                   family$varfun)       # quasi
    if (!(vTyp <- match(vNam, varNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    c(fTyp = fTyp, lTyp = lTyp, vTyp = vTyp)
}

convergenceMessage <- function(cvg)
### Create the convergence message
{
    msg <- switch(as.character(cvg),
                  "3" = "X-convergence (3)",
                  "4" = "relative convergence (4)",
                  "5" = "both X-convergence and relative convergence (5)",
                  "6" = "absolute function convergence (6)",

                  "7" = "singular convergence (7)",
                  "8" = "false convergence (8)",
                  "9" = "function evaluation limit reached without convergence (9)",
                  "10" = "iteration limit reached without convergence (9)",
                  "14" = "storage has been allocated (?) (14)",

                  "15" = "LIV too small (15)",
                  "16" = "LV too small (16)",
                  "63" = "fn cannot be computed at initial par (63)",
                  "65" = "gr cannot be computed at initial par (65)")
    if (is.null(msg))
        msg <- paste("See PORT documentation.  Code (", cvg, ")", sep = "")
    msg
}

mer_finalize <- function(ans)
{
    .Call(mer_optimize, ans)
    if (ans@dims[["cvg"]] > 6) warning(convergenceMessage(ans@dims[["cvg"]]))
    # Call(mer_update_ranef, ans) # blme edit: not needed, added to mer_optimize
    .Call(mer_update_mu, ans)
    ans
}
