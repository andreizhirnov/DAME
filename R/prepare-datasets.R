
makeframes <- function(data, f, bins, at=NULL, weights=NULL){
## extract all variables  
  allvars <- all.vars(f) ## all variables included in the formula 
  at_vars <- names(at)
  at_vars_n <- paste0("at.", at_vars)
  varying_vars <- setdiff(allvars,names(at))
  varying_vars_n <- paste0("va.",varying_vars) 
## fixed values
  if (length(at_vars)>0) {
    at_grid <- at_dt <- expand.grid(as.list(at)) 
    colnames(at_grid) <- at_vars_n
    at_grid[["at_id"]] <- seq_len(nrow(at_grid))
    if (length(varying_vars)==0 && length(at_vars) >0) {
      coord <- seq_len(nrow(at_dt))-1L
      wei <- data.frame(x = coord, y = coord, v = 1)
      at_grid[["at_id"]] <- NULL
      return(list(samples = at_dt, 
                  grid = at_grid, 
                  wei_locs = as.matrix(wei[,c("x","y")]) |> t(),
                  wei_vals = wei$v))
    }
  }
  
## dataset with varying information  
  usable <- stats::complete.cases(data[setdiff(allvars,names(at))])
  varying <- subset(data, 
                    subset=usable, 
                    select=varying_vars)
  colnames(varying) <- varying_vars_n
  varying[["bin"]] <- bins[["bin"]][usable]
  
## compressed
  key <- do.call(paste, c(varying, sep="\r"))
  varying$row_id <- as.integer(factor(key))
  compressed <- varying[!duplicated(varying$row_id),]
  compressed <- compressed[order(compressed$row_id),]
  
## wt
  wei <- data.frame(x = compressed$row_id-1L, 
                    y = compressed$bin-1L)
  if (is.null(weights)) { 
    wei[["v"]] <- tabulate(varying$row_id) |> as.vector()
  } else {
    wei[["v"]] <- rowsum(weights[usable], varying$row_id, reorder=FALSE)[,1] |> as.vector()
  }
  
# samples  
  samples <- subset(compressed, select = varying_vars_n)
  grid <- bins[["xwalk"]][match(bins[["xwalk"]][["bin"]], compressed$bin),]
  
# no at-variables
  if (length(at_vars)==0) {
    colnames(samples) <- substring(colnames(samples),4L)
    return(list(samples=samples, 
                grid=grid, 
                wei_locs = as.matrix(wei[,c("x","y")]) |> t(),
                wei_vals = wei$v))
  }

# combination of varying variables and at values
  samples <- samples |> merge(at_grid, by=NULL)
  samples[["at_id"]] <- NULL
  colnames(samples) <- substring(colnames(samples),4L)
  grid <- grid |> merge(at_grid, by=NULL)
  mx <- max(wei$x) + 1L
  my <- max(wei$y) + 1L 
  wei <- wei |> merge(at_grid[,"at_id",drop=FALSE], by=NULL) |> within({
    x <- mx*(at_id-1L)+x 
    y <- my*(at_id-1L)+y 
  })
  return(list(samples=samples, 
              grid=grid, 
              wei_locs = as.matrix(wei[,c("x","y")]) |> t(),
              wei_vals = wei$v))   
}
