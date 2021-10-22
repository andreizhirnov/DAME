#' @importFrom stats aggregate na.omit complete.cases quantile median

makeframes.dame <- function(data, allvars, at=NULL, bin_id, weights) {
  at_dt <- expand.grid(as.list(at))
  if (length(setdiff(allvars,names(at)))==0) {
    grid <- at_dt
    colnames(grid) <- paste0("at.", colnames(grid))
    data.compressed <- at_dt
    wmat <- diag(nrow(at_dt))
    return(list(data=data.compressed, grid=grid, wmat=wmat))
  }
  usable <- complete.cases(data[setdiff(allvars,names(at))])
  varying <- subset(data, usable, select=setdiff(allvars,names(at)))
  colnames(varying) <- paste0("va.",colnames(varying))
  varying[["bin_id"]] <- bin_id[usable]
  uni_row <- which(!duplicated(varying))
  uni_bin <- sort(unique(varying[["bin_id"]]))
  # single bin
  if (length(uni_bin)==1 && nrow(at_dt) <= 1) {
    grid <- data.frame(bin_id=uni_bin)
    data.compressed <- varying[uni_row,,drop=FALSE]
    data.compressed[["row_id"]] <- seq_along(uni_row)
    varying[["weight"]] <- weights[usable]
    temp <- merge(varying, data.compressed)
    counts <- aggregate(list("X_num" = temp$weight), by = temp[c("row_id","bin_id")], FUN = sum)
    wmat <- as.matrix(counts[["X_num"]]/sum(counts[["X_num"]]))
    if (ncol(at_dt)>0) {
      colnames(at_dt) <- paste0("at.", colnames(at_dt))
      data.compressed <- merge(data.compressed, at_dt, by=NULL)
      grid <- merge(grid, at_dt, by=NULL)
      }
    data.compressed[["row_id"]] <- data.compressed[["bin_id"]] <- NULL
    colnames(data.compressed) <- substring(colnames(data.compressed),4)
    return(list(data=data.compressed, grid=grid, wmat=wmat))
  }
  # everything else
  grid <- data.frame(bin_id=uni_bin, grid_id=seq_along(uni_bin))
  ##  grid$grid_id <- seq_len(nrow(grid))
  data.compressed <- varying[uni_row,]
  data.compressed[["row_id"]] <- seq_len(nrow(data.compressed))
  varying[["weight"]] <- weights[usable]
  temp <- merge(varying, data.compressed)
  counts <- aggregate(list("X_num" = temp$weight), by = temp[c("row_id","bin_id")], FUN = sum)
  if (ncol(at_dt)==0) {
    counts <- counts[order(counts$row_id),]
    counts$grid_id <- match(counts$bin_id, uni_bin)
  } else {
    colnames(at_dt) <- paste0("at.", colnames(at_dt))
    at_dt$at_id <- seq_len(nrow(at_dt))
    data.compressed <- merge(data.compressed, at_dt,by=NULL)
    data.compressed <- data.compressed[order(data.compressed$row_id, data.compressed$at_id),]
    grid <- merge(grid, at_dt, by=NULL)
    grid$grid_id <- paste0(grid$grid_id,":",grid$at_id)
    counts <- merge(counts, at_dt[,"at_id",drop=FALSE], by=NULL)
    counts <- counts[order(counts$row_id, counts$at_id),]
    counts$grid_id <- paste0(match(counts$bin_id, uni_bin),":",counts$at_id)
    data.compressed$at_id <- grid$at_id  <- NULL
  }
  # wmat <- matrix(nrow=nrow(counts), ncol=nrow(grid))
  # for (j in seq_len(nrow(grid))) {
  #   t <- counts$X_num*as.numeric(counts$grid_id==grid$grid_id[j])
  #   wmat[,j] <- t/sum(t)
  # }

  grid_diag <- diag(nrow(grid))
  colnames(grid_diag) <- paste0("X.",seq_len(nrow(grid)))
  combo <- merge(counts, data.frame(grid_id=grid$grid_id,grid_diag), by="grid_id")
  combo <- combo[order(combo$row_id),]
  grid_prep <- as.matrix(combo[["X_num"]] * combo[,paste0("X.",seq_len(nrow(grid))), drop=FALSE])
  wmat <- grid_prep %*% diag(1/colSums(grid_prep))

  grid$grid_id <- NULL
  data.compressed[["row_id"]] <- data.compressed[["bin_id"]] <- NULL
  colnames(data.compressed) <- substring(colnames(data.compressed),4)
  return(list(data=data.compressed, grid=grid, wmat=wmat))
}

makeframes.me <- function(data,tovary,at=NULL) {
  grid <- na.omit(data[,tovary,drop=FALSE])
  if (length(at) > 0) grid <- merge(expand.grid(as.list(at)), grid, by = NULL)
  grid
}

makeframes.mem <- function(at) {
  expand.grid(at)
}
