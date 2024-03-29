#' @importFrom methods as

make.dydm <- function (link)
{
  switch(EXPR=link,
         "logit" = {
           ym <- function(mu) 1/(1+exp(-mu))
           dydm <- function(mu) {d <- ym(mu); d*(1-d)}
           d2ydm2 <- function(mu) {d <- ym(mu); d*(1-d)*(1-2*d)}
         },
         "probit" = {
           ym <- function(mu) stats::pnorm(mu)
           dydm <- function(mu) stats::dnorm(mu)
           d2ydm2 <- function(mu) -mu*stats::dnorm(mu)
         },
         "cauchit" = {
           ym <- function(mu) stats::pcauchy(mu)
           dydm <- function(mu) stats::dcauchy(mu)
         },
         "cloglog" = {
           ym <- function(mu) -expm1(-exp(mu))
           dydm <- function(mu) exp(mu - exp(mu))
           d2ydm2 <- function(mu) -exp(mu - exp(mu))*expm1(mu)
         },
         "identity" = {
           ym <- function(mu) mu
           dydm <- function(mu) 1
           d2ydm2 <- function(mu) 0
         },
         "log" = {
           ym <- function(mu) exp(mu)
           dydm <- function(mu) exp(mu)
           d2ydm2 <- function(mu) exp(mu)
         },
         "sqrt" = {
           ym <- function(mu) mu^2
           dydm <- function(mu) 2*mu
           d2ydm2 <- function(mu) 2
         },
         "1/mu^2" = {
           ym <- function(mu) 1/sqrt(mu)
           dydm <- function(mu) -(mu^(-3/2))/2
           d2ydm2 <- function(mu) 3*(mu^(-5/2))/4
         },
         "inverse" = {
           ym <- function(mu) 1/mu
           dydm <- function(mu) -mu^(-2)
           d2ydm2 <- function(mu) 2*mu^(-3)
         }
  )
  list(ym=ym,dydm=dydm,d2ydm2=d2ydm2)
}

clean.calls <- function(x) {
  if (length(x) == 1) x
  else {
    r <- x
    if (r[[1]] ==  quote(I)) r <- r[[-1]]
    if (r[1] == quote(`:`())) r[1] <- quote(`*`())
    as.call(lapply(r, function(z) {if (length(z) == 1) {return(z)} else {clean.calls(z)}}))
  }
}

mx <- function(mmat, coefficients, offset) {
  if (length(offset) == nrow(mmat)) {
    sweep(coefficients %*% t(mmat),2,-offset)
  } else {
    coefficients %*% t(mmat)
  }
}

make.dmdx <- function(formula, bnames, xvarname) {
  mx.parts <- lapply(bnames, function(x) clean.calls(str2lang(noquote(x))))
  names(mx.parts) <- bnames
  dmdx.parts <- lapply(mx.parts, function(x) {
    tryCatch(stats::D(x,xvarname), error = function(e) {cat("Warning: Could not find the derivative of",x," and will replace it with zero"); 0})
  })
  return(list(dmdx.parts=dmdx.parts, mx.parts=mx.parts))
}

dmdb <- function(mmat, data, mx.parts) {
  mmatd <- as.data.frame(mmat)
  upddata <- c(list(Intercept=1),mmatd[setdiff(names(mmatd),names(data))], data)
  do.call("cbind",lapply(mx.parts, eval, envir = upddata))
}
d2mdxdb <- function(mmat, data, dmdx.parts) {
  mmatd <- as.data.frame(mmat)
  upddata <- c(list(Intercept=1),mmatd[setdiff(names(mmatd),names(data))], data)
  do.call("cbind",lapply(dmdx.parts, eval, envir = upddata))
}
dmdx <- function(mmat, data, coefficients, dmdx.parts) {
  mmatd <- as.data.frame(mmat)
  upddata <- c(list(Intercept=1),mmatd[setdiff(names(mmatd),names(data))], data)
  parts <- do.call("cbind",lapply(dmdx.parts, eval, envir = upddata))
  coefficients %*% t(parts)
}

find.central <- function(x,data,weights=NULL) {
  y <- data[[x]]
  if (is.numeric(y) && length(weights) != length(y)) {
    mean(y,na.rm=TRUE)
  } else if (is.numeric(y)) {
    sum(y*weights, na.rm=TRUE)/sum((1-is.na(y))*weights, na.rm=TRUE)
  } else if (length(weights) != length(y)) {
    y <- as.factor(y)
    lev <- levels(y)
    tab <- tabulate(match(y, lev))
    factor(lev[which.max(tab)], levels=lev)
  } else {
    y <- as.factor(y)
    lev <- levels(y)
    tab <- aggregate(weights ~ y, FUN=sum, na.action=NULL, na.rm=TRUE)
    factor(tab$y[which.max(tab$weights)], levels=lev)
  }
}

make.bins <- function(x, nbins) {
  qrs <- sort(unique(stats::quantile(x, seq(0, 1, by = 1/nbins), na.rm=TRUE)))
  cuts <- as.numeric(cut(x, qrs, include.lowest = TRUE))
  mp <- aggregate(x ~ cuts, FUN=stats::median, na.action=NULL, na.rm=TRUE)
  mp[match(cuts,mp$cuts) ,"x"]
}

simulated.me <- function(discrete, discrete_step=1, iter, coefficients, vcov, data, x, formula, ym, dydm, wmat = NULL, pct, ...) {
  mf <- stats::model.frame(formula = formula, data = data)
  offset <- stats::model.offset(mf)
  mmat <- stats::model.matrix(object = formula, data = mf)
  beta.names <- intersect(names(coefficients),colnames(mmat))
  if (length(beta.names) == 0) stop("Failed to link the supplied coefficients to the formula", call. = FALSE)
  mmat <- mmat[,beta.names,drop=FALSE]
  coef_vect <- matrix(coefficients[beta.names], nrow=1L)
  var_covar <- as.matrix(vcov)[beta.names,beta.names]
  coef_matrix <- MASS::mvrnorm(n = iter, mu = coef_vect, Sigma = var_covar, empirical = TRUE)
  if (discrete == TRUE) {
    data_offset <- data
    data_offset[[x]] <- data_offset[[x]] + discrete_step
    mf_offset <- stats::model.frame(formula = formula, data = data_offset)
    mmat_offset <- stats::model.matrix(object = formula, data = mf_offset)
    mmat_offset <- mmat_offset[,beta.names]
    bulk <- try(ym(mx(mmat = mmat_offset, coefficients = coef_matrix, offset = offset)) - ym(mx(mmat = mmat, coefficients = coef_matrix, offset = offset)), silent=TRUE)
    if (inherits(bulk,'try-error') & nrow(mmat) < iter) {
      bulk <- matrix(NA,nrow=iter,ncol=nrow(mmat))
      for (i in seq_len(nrow(mmat))) {
        bulk[,i] <- ym(mx(mmat = mmat_offset[i,,drop=FALSE], coefficients = coef_matrix, offset = offset)) - ym(mx(mmat = mmat[i,,drop=FALSE], coefficients = coef_matrix, offset = offset))
      }
    } else {
      bulk <- matrix(NA,nrow=iter,ncol=nrow(mmat))
      for (i in seq_len(iter)) {
        bulk[i,] <- ym(mx(mmat = mmat_offset, coefficients = coef_matrix[i,,drop=FALSE], offset = offset)) - ym(mx(mmat = mmat, coefficients = coef_matrix[i,,drop=FALSE], offset = offset))
      }
    }
  } else {
    dmli <- make.dmdx(formula = formula, bnames = beta.names, xvarname=x)
      bulk <- try(dydm(mx(mmat = mmat, coefficients = coef_matrix, offset = offset)) * dmdx(mmat = mmat, data=data, coefficients = coef_matrix, dmdx.parts=dmli$dmdx.parts), silent=TRUE)
      if (inherits(bulk,'try-error') & nrow(mmat) < iter) {
      bulk <- matrix(NA,nrow=iter,ncol=nrow(mmat))
      for (i in seq_len(nrow(mmat))) {
        bulk[,i] <- dydm(mx(mmat = mmat[i,,drop=FALSE], coefficients = coef_matrix, offset = offset)) * dmdx(mmat = mmat[i,,drop=FALSE], data=data[i,,drop=FALSE], coefficients = coef_matrix, dmdx.parts=dmli$dmdx.parts)
      }
      } else if (inherits(bulk,'try-error')) {
      bulk <- matrix(NA,nrow=iter,ncol=nrow(mmat))
      for (i in seq_len(iter)) {
        bulk[i,] <- dydm(mx(mmat = mmat, coefficients = coef_matrix[i,,drop=FALSE], offset = offset)) * dmdx(mmat = mmat, data=data, coefficients = coef_matrix[i,,drop=FALSE], dmdx.parts=dmli$dmdx.parts)
      }
    }
  }
  if (is.null(wmat)) {
    me_matrix <- bulk
  } else {
    me_matrix <- bulk %*% wmat
  }
  est <- colMeans(me_matrix)
  se <- apply(me_matrix, 2L, stats::sd)
  quantiles <- lapply(pct/100, function(x) apply(me_matrix, 2L, stats::quantile, probs = x, names = FALSE))
  data.frame(est=as.vector(est),se=se,quantiles)
}

analytical.me <- function(discrete, discrete_step=1, coefficients, vcov, data, x, formula, ym, dydm, d2ydm2, wmat = NULL, pct, ...) {
  mf <- stats::model.frame(formula = formula, data = data)
  offset <- stats::model.offset(mf)
  mmat <- stats::model.matrix(object = formula, data = mf)
  beta.names <- intersect(names(coefficients),colnames(mmat))
  if (length(beta.names) == 0) stop("Failed to link the supplied coefficients to the formula", call. = FALSE)
  mmat <- mmat[,beta.names,drop=FALSE]
  coef_vect <- matrix(coefficients[beta.names], nrow=1L)
  var_covar <- as.matrix(vcov)[beta.names,beta.names]
  dmli <- make.dmdx(formula = formula, bnames = beta.names, xvarname=x)
  m_0 <- as.vector(mx(mmat = mmat, coefficients = coef_vect, offset = offset))
  if (discrete == TRUE) {
    data_offset <- data
    data_offset[[x]] <- data_offset[[x]] + discrete_step
    mf_offset <- stats::model.frame(formula = formula, data = data_offset)
    mmat_offset <- stats::model.matrix(object = formula, data = mf_offset)
    mmat_offset <- mmat_offset[,beta.names,drop=FALSE]
    m_1 <- as.vector(mx(mmat = mmat_offset, coefficients = coef_vect, offset = offset))
    est <- ym(m_1) - ym(m_0)
    L <- dydm(m_1)*dmdb(mmat = mmat_offset, data=data_offset, mx.parts=dmli$mx.parts) - dydm(m_0)*dmdb(mmat = mmat, data=data, mx.parts=dmli$mx.parts)
    } else {
    est <- dydm(m_0)*as.vector(dmdx(mmat = mmat, data=data, coefficients = coef_vect, dmdx.parts=dmli$dmdx.parts))
    der2 <- d2mdxdb(mmat=mmat, data = data, dmdx.parts=dmli$dmdx.parts)
    if (nrow(der2)==1L && nrow(mmat) >1L) {
      der2 <- matrix(data = as.vector(der2), nrow = nrow(mmat), ncol = ncol(der2), byrow = TRUE)
      }
    L <- d2ydm2(m_0)*as.vector(dmdx(mmat = mmat, data=data, coefficients = coef_vect, dmdx.parts=dmli$dmdx.parts))*dmdb(mmat = mmat, data=data, mx.parts=dmli$mx.parts) + dydm(m_0)*der2
    }
  if (!is.null(wmat)) {
    est <- matrix(est,nrow=1L) %*% wmat
    L <- t(wmat) %*% L
  }
  var <- as.vector(apply(L,1L, function(x) (x %*% var_covar) %*% x))
  var[var < 0] <- 0
  se <- sqrt(var)
  quantiles <- lapply(pct/100, stats::qnorm, mean = as.vector(est), sd = se)
  data.frame(est=as.vector(est),se=se,quantiles)
}

check.required <- function(name, type, list = NULL) {
  if (is.null(list)) {
    if (!exists(name, envir=parent.frame(), inherits=FALSE)) {
      stop(paste0("Required argument '",name[1],"' is missing"), call. = FALSE)
    } else if (type=="character") {
      if (!inherits(eval(as.name(name), envir=parent.frame()),type)) {
        stop(paste0("Argument '",name[1],"' must be a character string"), call. = FALSE)
      }
    } else {
      candidate <- NULL
      candidate <- tryCatch(as(eval(as.name(name), envir=parent.frame()),type), error = function(e) return(NULL))
      if (is.null(candidate)) stop(paste0("Argument '",name[1],"' must be a ",type[1]), call. = FALSE)
    }
  } else {
    if (is.null(list[[name]])) {
      stop(paste0("Required argument '",name[1],"' is missing"), call. = FALSE)
    } else if (type=="character") {
      if (!inherits(list[[name]],type)) {
        stop(paste0("Argument '",name[1],"' must be a character string"), call. = FALSE)
      }
    } else {
      candidate <- NULL
      candidate <- tryCatch(as(list[[name]],type), error = function(e) return(NULL))
      if (is.null(candidate)) stop(paste0("Error: Argument '",name[1],"' must be a ",type[1]), call. = FALSE)
    }
  }
}
