make.dydm <- function (link)
{
  switch(EXPR=link,
         "logit" = {
           ym <- function(mu) 1/(1+exp(-mu))
           dydm <- function(mu) {d <- ym(mu); d*(1-d)}
           d2ydm2 <- function(mu) {d <- ym(mu); d*(1-d)*(1-2*d)}
         },
         "probit" = {
           ym <- function(mu) pnorm(mu)
           dydm <- function(mu) dnorm(mu)
           d2ydm2 <- function(mu) -mu*dnorm(mu)
         },
         "cauchit" = {
           ym <- function(mu) pcauchy(mu)
           dydm <- function(mu) dcauchy(mu)
           d2ydm2 <- function(mu) { d <- dcauchy(mu); -2*mu*d^2}
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
  assign("ym", ym, env=parent.frame())
  assign("dydm", dydm, env=parent.frame())
  assign("d2ydm2", d2ydm2, env=parent.frame())
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

mx <- function(mmat, coefficients) {
  coefficients %*% t(mmat)
}

make.dmdx <- function(formula, bnames, xvarname) {
  mx.parts <- lapply(bnames, function(x) clean.calls(str2lang(noquote(x))))
  names(mx.parts) <- bnames
  dmdx.parts <- lapply(mx.parts, function(x) {
    tryCatch(D(x,xvarname), error = function(e) {cat("Warning: Could not find the derivative of",x," and will replace it with zero"); 0})
  })
  dmdb <- function(mmat, data) {
    mmatd <- as.data.frame(mmat)
    upddata <- c(list(Intercept=1),mmatd[setdiff(names(mmatd),names(data))], data)
    do.call("cbind",lapply(mx.parts, function(x) eval(x, envir = upddata)))
  }
  d2mdxdb <- function(mmat, data) {
    mmatd <- as.data.frame(mmat)
    upddata <- c(list(Intercept=1),mmatd[setdiff(names(mmatd),names(data))], data)
    do.call("cbind",lapply(dmdx.parts, function(x) eval(x, envir = upddata)))
  }
  dmdx <- function(mmat, data, coefficients) {
    mmatd <- as.data.frame(mmat)
    upddata <- c(list(Intercept=1),mmatd[setdiff(names(mmatd),names(data))], data)
    parts <- do.call("cbind",lapply(dmdx.parts, function(x) eval(x, envir = upddata)))
    coefficients %*% t(parts)
  }
  assign("dmdx",dmdx, env=parent.frame())
  assign("dmdb", dmdb , env=parent.frame())
  assign("d2mdxdb", d2mdxdb, env=parent.frame())
}

find.mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

simulated.me <- function(discrete, discrete_step, iter, coefficients, variance, data, x, formula, ym, mx, dydm, wmat = NULL, pct) {
  mf <- model.frame(formula = formula, data = data)
  mmat <- model.matrix(object = formula, data = mf)
  beta.names <- intersect(names(coefficients),colnames(mmat))
  if (length(beta.names) == 0) stop("Failed to link the supplied coefficients to the formula", call. = FALSE)
  mmat <- mmat[,beta.names,drop=FALSE]
  coef_vect <- matrix(coefficients[beta.names], nrow=1L)
  var_covar <- variance[beta.names,beta.names]
  coef_matrix <- MASS::mvrnorm(n = iter, mu = coef_vect, Sigma = var_covar, empirical = TRUE)
  if (discrete == TRUE) {
    data_offset <- data
    data_offset[[x]] <- data_offset[[x]] + discrete_step
    mf_offset <- model.frame(formula = formula, data = data_offset)
    mmat_offset <- model.matrix(object = formula, data = mf_offset)
    mmat_offset <- mmat_offset[,beta.names]
    bulk <- try(ym(mx(mmat = mmat_offset, coefficients = coef_matrix)) - ym(mx(mmat = mmat, coefficients = coef_matrix)), silent=TRUE)
    if (inherits(bulk,'try-error') & nrow(mmat) < iter) {
      bulk <- matrix(NA,nrow=iter,ncol=nrow(mmat))
      for (i in 1L:nrow(mmat)) {
        bulk[,i] <- ym(mx(mmat = mmat_offset[i,,drop=FALSE], coefficients = coef_matrix)) - ym(mx(mmat = mmat[i,,drop=FALSE], coefficients = coef_matrix))
      }
    } else {
      bulk <- matrix(NA,nrow=iter,ncol=nrow(mmat))
      for (i in 1L:iter) {
        bulk[i,] <- ym(mx(mmat = mmat_offset, coefficients = coef_matrix[i,,drop=FALSE])) - ym(mx(mmat = mmat, coefficients = coef_matrix[i,,drop=FALSE]))
      }
    }
  } else {
    make.dmdx(formula = formula, bnames = beta.names, xvarname=x)
      bulk <- try(dydm(mx(mmat = mmat, coefficients = coef_matrix)) * dmdx(mmat = mmat, data=data, coefficients = coef_matrix), silent=TRUE)
      if (inherits(bulk,'try-error') & nrow(mmat) < iter) {
      bulk <- matrix(NA,nrow=iter,ncol=nrow(mmat))
      for (i in 1L:nrow(mmat)) {
        bulk[,i] <- dydm(mx(mmat = mmat[i,,drop=FALSE], coefficients = coef_matrix)) * dmdx(mmat = mmat[i,,drop=FALSE], data=data[i,,drop=FALSE], coefficients = coef_matrix)
      }
      } else if (inherits(bulk,'try-error')) {
      bulk <- matrix(NA,nrow=iter,ncol=nrow(mmat))
      for (i in 1L:iter) {
        bulk[i,] <- dydm(mx(mmat = mmat, coefficients = coef_matrix[i,,drop=FALSE])) * dmdx(mmat = mmat, data=data, coefficients = coef_matrix[i,,drop=FALSE])
      }
    }
  }
  if (is.null(wmat)) {
    me_matrix <- bulk
  } else {
    me_matrix <- bulk %*% wmat
  }
  est <- colMeans(me_matrix)
  se <- apply(me_matrix, 2L, sd)
  quantiles <- lapply(pct/100, function(x) apply(me_matrix, 2L, quantile, probs = x, names = FALSE))
  data.frame(est=as.vector(est),se=se,quantiles)
}

analytical.me <- function(discrete, discrete_step, coefficients, variance, data, x, formula, ym, mx, dydm, d2ydm2, wmat = NULL, pct) {
  mf <- model.frame(formula = formula, data = data)
  mmat <- model.matrix(object = formula, data = mf)
  beta.names <- intersect(names(coefficients),colnames(mmat))
  if (length(beta.names) == 0) stop("Failed to link the supplied coefficients to the formula", call. = FALSE)
  mmat <- mmat[,beta.names,drop=FALSE]
  coef_vect <- matrix(coefficients[beta.names], nrow=1L)
  var_covar <- variance[beta.names,beta.names]
  make.dmdx(formula = formula, bnames = beta.names, xvarname=x)
  m_0 <- as.vector(mx(mmat = mmat, coefficients = coef_vect))
  if (discrete == TRUE) {
    data_offset <- data
    data_offset[[x]] <- data_offset[[x]] + discrete_step
    mf_offset <- model.frame(formula = formula, data = data_offset)
    mmat_offset <- model.matrix(object = formula, data = mf_offset)
    mmat_offset <- mmat_offset[,beta.names,drop=FALSE]
    m_1 <- as.vector(mx(mmat = mmat_offset, coefficients = coef_vect))
    est <- ym(m_1) - ym(m_0)
    L <- dydm(m_1)*dmdb(mmat = mmat_offset, data=data_offset) - dydm(m_0)*dmdb(mmat = mmat, data=data)
    } else {
    est <- dydm(m_0)*as.vector(dmdx(mmat = mmat, data=data, coefficients = coef_vect))
    L <- d2ydm2(m_0)*as.vector(dmdx(mmat = mmat, data=data, coefficients = coef_vect))*dmdb(mmat = mmat, data=data) + dydm(m_0)*d2mdxdb(mmat=mmat, data = data)
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

check.args <- function(args,checks) {
  errors.required <- setdiff(checks[["required"]],names(args))
  if (length(errors.required)>0) {
    stop(paste0("Required argument '",errors.required[1],"' is missing"), call. = FALSE)
  }
  tocheck <- intersect(names(checks[["types"]]),names(args))
  if (length(tocheck)>0) {
    for (i in tocheck){
      candidate <- NULL
      candidate <- tryCatch(as(eval(args[[i]]),checks[["types"]][[i]]),
                            error = function(e) return(ifelse(checks[["types"]][[i]] == "character", toString(args[[i]]), NULL)))
      if (is.null(candidate)) stop(paste0("Error: Argument '",i,"' is not ",checks[["types"]][[i]]), call. = FALSE)
      if (is.null(checks[["lengths"]][[i]])) {
        assign(i, candidate, env=parent.frame())
      } else if (length(candidate) == checks[["lengths"]][[i]]) {
        assign(i,candidate, env=parent.frame())
      } else {
        stop(paste0("Error: Argument '",i,"' must be of length ",checks[["lengths"]][[i]]), call. = FALSE)
      }
    }
  }
}
