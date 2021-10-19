#' @importFrom stats aggregate na.omit complete.cases quantile median

makeframes.dame <- function(data, allvars, at=NULL, bin_id, weight) {
  at_dt <- expand.grid(as.list(at))
  usable <- complete.cases(data[setdiff(allvars,names(at))])
  varying <- subset(data, usable, select=setdiff(allvars,names(at)))
  colnames(varying) <- paste0("va.",colnames(varying))
  varying[["bin_id"]] <- bin_id[usable]
  varying[["weight"]] <- weight[usable]
  uni_row <- which(!duplicated(varying))
  uni_bin <- sort(unique(varying[["bin_id"]]))

# single bin
  if (length(uni_bin)==1 && nrow(at_dt) <= 1) {
    grid <- data.frame(bin_id=uni_bin)
    data.compressed <- varying[uni_row,]
    data.compressed[["row_id"]] <- seq_along(uni_row)
    temp <- merge(varying, data.compressed)
    counts <- aggregate(list("X_num" = temp$weight), by = temp[c("row_id","bin_id")], FUN = sum)
    wmat <- as.matrix(counts[["X_num"]]/sum(counts[["X_num"]]))
    data.compressed[["row_id"]] <- data.compressed[["bin_id"]] <- data.compressed[["at_id"]] <- NULL
    colnames(data.compressed) <- substring(colnames(data.compressed),4)
    return(list(data.compressed=data.compressed, grid=grid, wmat=wmat))
  }
# everything else
  grid <- data.frame(bin_id=uni_bin, grid_id=seq_along(uni_bin))
##  grid$grid_id <- seq_len(nrow(grid))
  data.compressed <- varying[uni_row,]
  data.compressed[["row_id"]] <- seq_len(nrow(data.compressed))
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
    grid$grid_id <- paste0(grid_id,":",at_id)
    counts <- merge(counts, at_dt[,"at_id",drop=FALSE], by=NULL)
    counts <- counts[order(counts$row_id, counts$at_id),]
    counts$grid_id <- paste0(counts$grid_id,":",counts$at_id)
 ##   counts <- merge(counts, grid[,c("bin_id","at_id","grid_id"),drop=FALSE], by=c("bin_id","at_id"))
    data.compressed$at_id <- grid$at_id  <- NULL
  }
  # wmat <- matrix(nrow=nrow(counts), ncol=nrow(grid))
  # for (j in grid$grid_id) {
  #   t <- counts$X_num*as.numeric(counts$grid_id==j)
  #   wmat[,j] <- t/sum(t)
  # }

  grid_diag <- diag(nrow(grid))
  colnames(grid_diag) <- paste0("X.",seq_len(nrow(grid)))
  combo <- merge(counts, cbind(grid_id=grid$grid_id,grid_diag), by="grid_id")
  combo <- combo[order(combo$row_id),]
  grid_prep <- as.matrix(combo[["X_num"]] * combo[,paste0("X.",seq_len(nrow(grid)))])
  wmat <- grid_prep %*% diag(1/colSums(grid_prep))

  grid$grid_id <- NULL
  data.compressed[["row_id"]] <- data.compressed[["bin_id"]] <- NULL
  colnames(data.compressed) <- substring(colnames(data.compressed),4)
  return(list(data.compressed=data.compressed, grid=grid, wmat=wmat))
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
