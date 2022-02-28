#' Max-Ratio Projection function
#' @description function to project high-dimensional datasets to lower dimention with max-ratio projection.
#' @param data The dataset to apply MRP. Each row is an observation.
#' @param cl The class identification for each observation. The length of \code{cl} should be the same as the number of rows of \code{data}.
#' @param nproj The number of max-ratio directions to be used in projecting the original data to the projected data.
#' @param message Logical. Wheather to show the accumulative variance explained by the projection directions or not.
#' @return A list with the elements
#' \item{projected_df}{The projected data with selected number of max-ratio directions.}
#' \item{pccumvar}{The cummulative variance explained by the max-ratio principal components.}
#' @export
mrp <- function(data, cl, nproj = 4, message = T, ...){
  n <- nrow(data)
  p <- ncol(data)
  m <- nlevels(cl)
  d <- min(p, min(summary(cl)))
  
  if (p < nproj) {
    stop(paste0("Asking for ", nproj, " directions, but only ", p, " dimensions."))
  }
  
  if (d < nproj) {
    stop(paste0("Asking for ", nproj, " directions, but minimum of ", d, " observation(s) in a cluster."))
  }
  
  class <- as.factor(as.numeric(cl))
  df <- data.frame(data,class)
  
  # Step 1  
  rotation <- matrix(0,ncol = d,nrow = p)
  for(i in 1:m){
    g.dat <- df[df$class == i,-(p+1)]
    g.dat.c <- apply(g.dat,2,function(x) x-mean(x))
    
    if(d > p) rotation <- rotation + eigen(cov(g.dat))$vectors
    else rotation <- rotation + svd(t(g.dat))$u[,1:d]
  }
  
  df.svd <- svd(rotation)
  df.rotation <- df.svd$u %*% t(df.svd$v)
  dat <- as.matrix(df[,-(p+1)]) %*% df.rotation
  
  # Step 2
  SST <- cov(dat)*(n - 1)
  df.t <- data.frame(dat,class) 
  SSG <- matrix(0,ncol = d,nrow = d)
  for(i in 1:m){
    g.dat <- df.t[df.t$class == i,-(d+1)]
    SSG <- SSG + cov(g.dat)*(nrow(g.dat) - 1)
  }
  SSB <- SST - SSG 
  SST <- (SST + t(SST))/2
  eig <- eigen(SST)
  d <- eig$values
  dd_ <- sqrt(1/d)
  dd_[is.na(dd_)] <- 0
  SST.sqrtinv <- eig$vectors %*% diag(sqrt(dd_)) %*% t(eig$vectors)
  W <- SST.sqrtinv%*%SSB%*%SST.sqrtinv
  W <- (W+t(W))/2
  r <- svd(W)
  
  
  pccumvar <- cumsum((r$d))/sum((r$d))
  if (message){
    cat("cumulative variance explained:", pccumvar, "\n")
  }
  if (nproj > length(r$d)){
    stop("Cannot get more max-ratio directions than ", length(r$d))
  }
  w <- r$v[,1:nproj]
  rotation1 <- SST.sqrtinv%*%w
  dat.rotation <- sweep(rotation1, MARGIN = 2, STATS = sqrt(colSums(rotation1^2)), FUN = "/")
  dat.red <- as.matrix(dat) %*% dat.rotation
  df1 <- data.frame(dat.red)
  colnames(df1) <- paste0("M[", 1:length(df1), "]")
  
  return(list(projected_df = df1, pccumvar = pccumvar))
}
