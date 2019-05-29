# Functions to manipulate the formula and retrieve random effect information


safeDeparse <- function(barsObj, collapse=" ") paste(deparse(barsObj, 500L), collapse=collapse)
# safeDeparse = direct copy from lme4


barnames <- function(bars) vapply(bars, function(x) safeDeparse(x[[3]]), "")
# barnames = direct copy from lme4


RHSForm <- function(form,as.form=FALSE) {
  # RHSForm = direct copy from lme4
  rhsf <- form[[length(form)]]
  if (as.form) reformulate(deparse(rhsf)) else rhsf
}


makeFac <- function(x,char.only=FALSE) {
  # makeFac = direct copy from lme4
  if (!is.factor(x) && (!char.only || is.character(x))) factor(x) else x
}

factorize <- function(x,frloc,char.only=FALSE) {
  # factorize = direct copy from lme4
  ## convert grouping variables to factors as necessary
  for (i in all.vars(RHSForm(x))) {
    if (!is.null(curf <- frloc[[i]]))
      frloc[[i]] <- makeFac(curf,char.only)
  }
  return(frloc)
}

##' @importFrom Matrix sparseMatrix
##' @importMethodsFrom Matrix t diag
##' @export
mkBlist_glmmPen <- function(x,frloc, drop.unused.levels=TRUE,
                            reorder.vars=FALSE) {
  # Slightly modified version of mkBlist from lme4
  ## Changes from original: sm <- KhatriRao(...) line 
  frloc <- factorize(x,frloc)
  ## try to evaluate grouping factor within model frame ...
  if (is.null(ff <- tryCatch(eval(substitute(makeFac(fac),
                                             list(fac = x[[3]])), frloc),
                             error=function(e) NULL)))
    stop("couldn't evaluate grouping factor ",
         deparse(x[[3]])," within model frame:",
         " try adding grouping factor to data ",
         "frame explicitly if possible",call.=FALSE)
  if (all(is.na(ff)))
    stop("Invalid grouping factor specification, ",
         deparse(x[[3]]),call.=FALSE)
  ## NB: *also* silently drops <NA> levels - and mkReTrms() and hence
  ##     predict.merMod() have relied on that property  :
  if (drop.unused.levels) ff <- factor(ff, exclude=NA)
  nl <- length(levels(ff))
  ## this section implements eq. 6 of the JSS lmer paper
  ## model matrix based on LHS of random effect term (X_i)
  ##    x[[2]] is the LHS (terms) of the a|b formula
  mm <- model.matrix(eval(substitute( ~ foo, list(foo = x[[2]]))), frloc)
  if (reorder.vars) {
    mm <- mm[colSort(colnames(mm)),]
  }
  ## this is J^T (see p. 9 of JSS lmer paper)
  ## construct indicator matrix for groups by observations
  ## use fac2sparse() rather than as() to allow *not* dropping
  ## unused levels where desired
  sm <- fac2sparse(ff, to = "d",
                   drop.unused.levels = drop.unused.levels)
  ### Change from original: switched order of matrices in KhariRao function
  sm <- KhatriRao(t(mm), sm)
  dimnames(sm) <- list(
    rep(levels(ff),each=ncol(mm)),
    rownames(mm))
  list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
}

##' @importFrom Matrix sparseMatrix drop0
##' @importMethodsFrom Matrix coerce
##' @export
mkReTrms_glmmPen <- function(bars, fr, drop.unused.levels=TRUE,
                             reorder.terms=TRUE,
                             reorder.vars=FALSE) {
  # Slightly modified version of mkReTrms from lme4
  ## Changes from original: call mkBlist_glmmPen instead of calling mkBlist
  if (!length(bars))
    stop("No random effects terms specified in formula",call.=FALSE)
  stopifnot(is.list(bars), vapply(bars, is.language, NA),
            inherits(fr, "data.frame"))
  names(bars) <- barnames(bars)
  term.names <- vapply(bars, safeDeparse, "")
  ## get component blocks
  ### Change from original mkReTrms: call mkBlist_glmmPen instead of call mkBlist
  blist <- lapply(bars, mkBlist_glmmPen, fr, drop.unused.levels,
                  reorder.vars = reorder.vars)
  nl <- vapply(blist, `[[`, 0L, "nl")   # no. of levels per term
  # (in lmer jss:  \ell_i)
  
  ## order terms stably by decreasing number of levels in the factor
  if (reorder.terms) {
    if (any(diff(nl) > 0)) {
      ord <- rev(order(nl))
      blist      <- blist     [ord]
      nl         <- nl        [ord]
      term.names <- term.names[ord]
    }
  }
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(rbind, Ztlist)  ## eq. 7, JSS lmer paper
  names(Ztlist) <- term.names
  q <- nrow(Zt)
  
  ## Create and install Lambdat, Lind, etc.  This must be done after
  ## any potential reordering of the terms.
  cnms <- lapply(blist, `[[`, "cnms")   # list of column names of the
  # model matrix per term
  nc <- lengths(cnms)                   # no. of columns per term
  # (in lmer jss:  p_i)
  nth <- as.integer((nc * (nc+1))/2)    # no. of parameters per term
  # (in lmer jss:  ??)
  nb <- nc * nl                         # no. of random effects per term
  # (in lmer jss:  q_i)
  ## eq. 5, JSS lmer paper
  if (sum(nb) != q) {
    stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)",
                 sum(nb),q))
  }
  boff <- cumsum(c(0L, nb))             # offsets into b
  thoff <- cumsum(c(0L, nth))           # offsets into theta
  
  Lambdat <-
    Matrix::t(do.call(sparseMatrix,
              do.call(rbind,
                      lapply(seq_along(blist), function(i)
                      {
                        mm <- matrix(seq_len(nb[i]), ncol = nc[i],
                                     byrow = TRUE)
                        dd <- diag(nc[i])
                        ltri <- lower.tri(dd, diag = TRUE)
                        ii <- row(dd)[ltri]
                        jj <- col(dd)[ltri]
                        ## unused: dd[cbind(ii, jj)] <- seq_along(ii)
                        data.frame(i = as.vector(mm[, ii]) + boff[i],
                                   j = as.vector(mm[, jj]) + boff[i],
                                   x = as.double(rep.int(seq_along(ii),
                                                         rep.int(nl[i], length(ii))) +
                                                   thoff[i]))
                      }))))
  thet <- numeric(sum(nth))
  ll <- list(Zt = drop0(Zt), theta = thet, Lind = as.integer(Lambdat@x),
             Gp = unname(c(0L, cumsum(nb))))
  ## lower bounds on theta elements are 0 if on diagonal, else -Inf
  ll$lower <- -Inf * (thet + 1)
  ll$lower[unique(diag(Lambdat))] <- 0
  ll$theta[] <- is.finite(ll$lower) # initial values of theta are 0 off-diagonal, 1 on
  Lambdat@x[] <- ll$theta[ll$Lind]  # initialize elements of Lambdat
  ll$Lambdat <- Lambdat
  # massage the factor list
  fl <- lapply(blist, `[[`, "ff")
  # check for repeated factors
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else asgn <- seq_along(fl)
  names(fl) <- ufn
  ## DON'T need fl to be a data.frame ...
  ## fl <- do.call(data.frame, c(fl, check.names = FALSE))
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist
  ll
} 
