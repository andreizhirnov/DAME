#' @importFrom stats aggregate na.omit


makeframes.ame <- function(data,allvars,at) {
  varying <- na.omit(data[intersect(allvars,setdiff(names(data),names(at)))])
  grid <- expand.grid(at)
  if (nrow(grid)== 0) {
    data.compressed <- aggregate(list("X_numerator" = rep(1/nrow(varying),nrow(varying))), by = varying, FUN = sum)
    wmat <- as.matrix(data.compressed[["X_numerator"]])
  } else if (nrow(grid) == 1) {
    data.compressed <- data.frame(at,aggregate(list("X_numerator" = rep(1/nrow(varying),nrow(varying))), by = varying, FUN = sum))
    wmat <- as.matrix(data.compressed[["X_numerator"]])
  } else {
    varying.compressed <- aggregate(list("X_numerator" = 1L:nrow(varying)), by = varying, FUN = length)
    data.compressed <- merge(expand.grid(at), varying.compressed, by = NULL)
    data.compressed[["X_rows"]] <- 1:nrow(data.compressed)
    grid.diag <- diag(nrow(grid))
    colnames(grid.diag) <- paste0("X_cols.",1L:ncol(grid.diag))
    combo <- merge(data.compressed, cbind(grid,grid.diag), by = names(grid))
    grid.prep <- as.matrix(combo[order(combo[["X_rows"]]),paste0("X_cols.",1L:nrow(grid))]) * combo[order(combo[["X_rows"]]),"X_numerator"]
    wmat <- grid.prep %*% diag(1/colSums(grid.prep))
  }
  list(data.compressed=data.compressed, grid=grid, wmat=wmat)
}

makeframes.dame <- function(data,allvars,at,over,x) {
  atmeans <- lapply(data[intersect(allvars,setdiff(names(data),c(names(at),x,over)))], function(x) if (is.numeric(x)) return(mean(x,na.rm=TRUE)) else return(find.mode(x)))
  varying <- na.omit(data[c(x,over)])
  if (length(over)==0) {
    grid <- expand.grid(c(at,atmeans))
  } else if (length(c(at,atmeans))==0) {
    grid <- unique(varying[over])
  } else {
    uniqover <- unique(varying[over])
    grid <- merge(expand.grid(c(at,atmeans)), uniqover, by = NULL)
  }
  if (nrow(grid)== 0) {
    data.compressed <- aggregate(list("X_numerator" = rep(1/nrow(varying),nrow(varying))), by = varying, FUN = sum)
    wmat <- as.matrix(data.compressed[["X_numerator"]])
  } else if (nrow(grid)== 1) {
    varying.compressed <- aggregate(list("X_numerator" = 1L:nrow(varying)), by = varying, FUN = length)
    data.compressed <- merge(grid, varying.compressed, by = over)
    data.compressed[["X_rows"]] <- 1:nrow(data.compressed)
    grid.diag <- diag(nrow(grid))
    colnames(grid.diag) <- paste0("X_cols.",1L:ncol(grid.diag))
    combo <- merge(data.compressed, cbind(grid,grid.diag), by = names(grid))
    grid.prep <- as.matrix(combo[order(combo[["X_rows"]]),paste0("X_cols.",1L:nrow(grid))]) * combo[order(combo[["X_rows"]]),"X_numerator"]
    wmat <- grid.prep / sum(grid.prep)
  } else {
    varying.compressed <- aggregate(list("X_numerator" = 1L:nrow(varying)), by = varying, FUN = length)
    data.compressed <- merge(grid, varying.compressed, by = over)
    data.compressed[["X_rows"]] <- 1:nrow(data.compressed)
    grid.diag <- diag(nrow(grid))
    colnames(grid.diag) <- paste0("X_cols.",1L:ncol(grid.diag))
    combo <- merge(data.compressed, cbind(grid,grid.diag), by = names(grid))
    grid.prep <- as.matrix(combo[order(combo[["X_rows"]]),paste0("X_cols.",1L:nrow(grid))]) * combo[order(combo[["X_rows"]]),"X_numerator"]
    wmat <- grid.prep %*% diag(1/colSums(grid.prep))
  }
  list(data.compressed=data.compressed, grid=grid, wmat=wmat)
}

makeframes.me <- function(data,allvars,at,over,x) {
  atmeans <- lapply(data[intersect(allvars,setdiff(names(data),c(names(at),x,over)))], function(x) if (is.numeric(x)) return(mean(x,na.rm=TRUE)) else return(find.mode(x)))
  varying <- na.omit(data[c(x,over)])
  if (length(c(at,atmeans)) == 0) {
    grid <- unique(varying[c(x,over)])
  } else {
    uniqover <- unique(varying[c(x,over)])
    grid <- merge(expand.grid(c(at,atmeans)), uniqover, by = NULL)
  }
  list(data.compressed=grid, grid=grid, wmat=NULL)
}

makeframes.mem <- function(data,allvars,at) {
  atmeans <- lapply(data[intersect(allvars,setdiff(names(data),names(at)))], function(x) if (is.numeric(x)) return(mean(x,na.rm=TRUE)) else return(find.mode(x)))
  grid <- expand.grid(c(at,atmeans))
  list(data.compressed=grid, grid=grid, wmat=NULL)
}
