#' @importFrom methods as

check_required <- function(name, type, list = NULL) {
  if (is.null(list)) {
    if (!exists(name, envir=parent.frame(), inherits=FALSE)) {
      stop(paste0("Required argument '",name[1L],"' is missing"), call. = FALSE)
    } else if (type=="character") {
      if (!inherits(eval(as.name(name), envir=parent.frame()),type)) {
        stop(paste0("Argument '",name[1L],"' must be a character string"), call. = FALSE)
      }
    } else {
      candidate <- NULL
      candidate <- tryCatch(as(eval(as.name(name), envir=parent.frame()),type), error = function(e) return(NULL))
      if (is.null(candidate)) stop(paste0("Argument '",name[1L],"' must be a ",type[1L]), call. = FALSE)
    }
  } else {
    if (is.null(list[[name]])) {
      stop(paste0("Required argument '",name[1L],"' is missing"), call. = FALSE)
    } else if (type=="character") {
      if (!inherits(list[[name]],type)) {
        stop(paste0("Argument '",name[1L],"' must be a character string"), call. = FALSE)
      }
    } else {
      candidate <- NULL
      candidate <- tryCatch(as(list[[name]],type), error = function(e) return(NULL))
      if (is.null(candidate)) stop(paste0("Error: Argument '",name[1L],"' must be a ",type[1L]), call. = FALSE)
    }
  }
}

check_link <- function(link_txt = NULL, model = NULL) {
  link <- link_txt[1L]
  if (is.null(link)) link <- eval(model)[["family"]][["link"]]
  check_required("link","character")
  
  l <- c('logit'=1L, 
         'probit'=2L, 
         'cauchit'=3L, 
         'cloglog'=4L, 
         'log'=5L,
         'sqrt'=6L, 
         "1/mu^2"=7L, 
         "inverse"=8L,
         "identity"=0L)[link]
  
  if (is.na(l)) {
    warning("Invalid link name. Valid links include 'logit','probit','cauchit',
         'cloglog','identity','log','sqrt','1/mu^2','inverse'. Defaulting to
            linear prediction.")
    return(0L)
  } else {
    return(l)
  }
}

check_data <- function(data=NULL, model=NULL) {
  d <- eval(data)
  if (is.null(d)) d <- eval(model)[["data"]]
  check_required("d","data.frame")
  return(d)
}

check_formula <- function(formula=NULL, model=NULL) {
  f <- formula
  if (is.null(f)) f <- stats::formula(model)
  f[[2L]] <- NULL
  check_required("f","formula")
  return(f)
}

check_weights <- function(weights=NULL, data=NULL) { 
  wt <- eval(weights)
  if (!is.null(wt) && !is.numeric(wt)) stop("'weights' must be a numeric vector", call. = FALSE)
  if (!is.null(wt) && length(wt) != nrow(wt)) stop("'weights' must have the same length as the dataset", call. = FALSE)
  return(wt)
}

check_at <- function(at = NULL, data){
  if (length(at)>0) {
    at <- as.list(at)
    for (v in names(at)) {
      if (is.character(at[[v]]) && !is.factor(at[[v]])) {
        xle <-  model[["xlevels"]][[v]]
        if (is.null(xle)) xle <- sort(unique(data[[v]]))
        if (is.null(xle)) {
          stop("Please convert the character variables in the 'at' list into factors", call. = FALSE)
        }
        if (any(!at[[v]] %in% xle)) {
          stop(paste0("Could not find all listed values of ",v," in the model"), call. = FALSE)
        }
        at[[v]] <- factor(at[[v]], levels=xle)
      }
    }
  }
  return(at)
}
 
clean_calls <- function(x) {
  if (length(x) == 1) x
  else {
    r <- x
    if (r[[1]] ==  quote(I)) r <- r[[-1]]
    if (r[1] == quote(`:`())) r[1] <- quote(`*`())
    as.call(lapply(r, function(z) {if (length(z) == 1) {return(z)} else {clean_calls(z)}}))
  }
}

find_central <- function(x, data, weights=NULL) {
  d <- data.frame(y = data[[x]])
  if (length(weights) != nrow(d)) {
    d$w <- 1 
  } else {
    d$w<- weights
  }
  d <- na.omit(d)
  if (is.numeric(d$y)) {
    return(with(d, sum(y*w)/sum(w)))
  } else {
    d$y <- as.factor(d$y)
    lev <- levels(d$y)
    return(with(d, tapply(w, y, FUN=sum)) |> 
             which.max() |> 
             names() |> 
             factor(levels=lev)) 
  }
}

find_central_by <- function(x, data, by, weights=NULL) {
  d <- data.frame(y = data[[x]], by = by)
  if (length(weights) != nrow(d)) {
    d$w <- 1 
  } else {
    d$w<- weights
  }
  d <- na.omit(d)
  if (is.numeric(d$y)) {
    agg <- with(d, tapply(y*w,by,FUN=sum) / tapply(w,by,FUN=sum))
    agg <- data.frame(by = names(agg), y=agg)
    colnames(agg) <- c(by, paste0("center.", x))
    return(agg)
  } else {
    d$y <- as.factor(d$y)
    lev <- levels(d$y)
    agg <- split(d, d$by) |>
      sapply(function(u) with(u, tapply(w, y, FUN=sum)) |> 
               which.max() |>  
               names()) |>
               factor(levels=lev) 
    agg <- data.frame(by = names(agg), y=agg)
    colnames(agg) <- c(by, paste0("center.", x))
    return(agg)
  }
}

make_bins <- function(bin_id=NULL, over=NULL, data=NULL, use_distinct_values = TRUE, nbins=NULL) {
  if (!is.null(bin_id)) {
    if (!is.numeric(bin_id)) stop("'bin_id' must be a numeric vector", call. = FALSE)
    if (length(bin_id) != nrow(data)) stop("'bin_id' must have the same length as the dataset", call. = FALSE)
    if (is.null(over)) {
      return(list(
          bin = as.integer(as.factor(bin_id)),
          xwalk = data.frame(bin_id = bin_id, 
                             bin = as.integer(as.factor(bin_id))) |> unique())
        )
    } else {
      t <- data.frame(bin_id = bin_id, bin = as.integer(as.factor(bin_id)))
      return(list(
        bin = t$bin,
        xwalk =  t |> unique() |>
          merge(find_central_by(x=over, data=data, by=bin_id, weights=weights), by="bin_id")
        ))
    }
  } else if (inherits(over, "character") && !is.null(data[[over]])) {
    if (use_distinct_values) { 
      t <- data[,over,drop=FALSE]
      colnames(t) <- paste0("center.", colnames(t))
      t$bin <- as.integer(as.factor(data[[over]]))
      return(list(
        bin = t$bin,
        xwalk = unique(t)
      ))
    } else {
      t <- data[,over,drop=FALSE]
      qrs <- sort(unique(stats::quantile(t[[1L]], seq(0, 1, by = 1/nbins), na.rm=TRUE)))
      t <- within(t, {
        bin <- as.integer(cut(t[[1L]], qrs, include.lowest = TRUE))
      })
      t[[paste0("center.", colnames(t))]] <- find_central_by(over, t, t$bin, weights)
      t[[over]] <- NULL
      return(list(
        bin = t$bin,
        xwalk = unique(t)
      ))
    }
    } else {
      return(list(
        bin = rep(1L, nrow(data)),
        xwalk = data.frame(bin=1L, bin_id=1L)
      )) 
    }
  }

make_mmat <- function(f, data) {
  mf <- stats::model.frame(formula = f, data = data) 
  stats::model.matrix(object = f, data = mf)
}

make_offset <- function(f, data) {
  stats::model.frame(formula = f, data = data) |> 
    stats::model.offset()
}

make_mmat_p <- function(f, data, x=NULL, discrete_step=NULL) {
  check_required("x","character")
  check_required(x, "numeric", list=data)
  check_required("discrete_step","numeric") 
  data[[x]] <- data[[x]] + discrete_step
  make_mmat(f, data)
}

make_d2mdxdb <- function(f, data, x=NULL) {
  check_required("x","character")
  check_required(x, "numeric", list=data)
  bnames <- stats::model.matrix(f, dt[0L,]) |> colnames()
  parts <- lapply(bnames, function(x) clean_calls(str2lang(noquote(x)))) 
  dmdx_parts <- lapply(setNames(parts, bnames), function(u) {
    tryCatch(stats::D(u,x), error = function(e) {cat("Warning: Could not find the derivative of",x," and will replace it with zero"); 0})
  }) 
  mmatd  <- make_mmat(f, data) |> as.data.frame()
  upddata <- c(list(Intercept=1),mmatd[setdiff(names(mmatd),names(data))], data)
  do.call("cbind",lapply(dmdx_parts, eval, envir = upddata))
}

check_coefs <- function(bnames, coefficients = NULL, model=NULL) {
  if (is.null(coefficients)) coefficients <- stats::coef(model)
  coefficients <- coefficients[bnames]
  check_required("coefficients", "numeric")
  return(coefficients)
}

check_vcov <- function(bnames, vcov = NULL, model=NULL) { 
if (is.null(vcov)) vcov <- stats::vcov(model)
vcov <- vcov[bnames, bnames]
check_required("vcov", "matrix")
return(vcov)
}

make_bounds <- function(pct) {
  check_required("pct", "numeric")
  if (any(pct > 100) || any(pct <0)) stop("Error: 'pct' must be between 0 and 100", call. = FALSE)
  p <- pct/100
  if (is.null(names(pct))) {
    names(p) <- paste0("p",pct)
  } else {
    names(p) <- make.names(names(pct))
  }
  return(p)
}
