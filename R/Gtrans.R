Fhat <- function(data){
  # returns a list of emperical cdfs. The list is of the same length of ncol(data)
  ecdfs <- list()
  
  for (j in 1:ncol(data)){
    ecdfs[[j]] <- ecdf(data[,j])
  }
  return(ecdfs)
}

epmf <- function(data){
  data;
  epm <- function(x) EnvStats::demp(x, data, discrete = T)
  return(epm)
}

pmf_hat <- function(data){
  # returns a list of emperical pmf functions
  epmfs <- list()
  for (j in 1:ncol(data)){
    #browser()
    epmfs <- c(epmfs, epmf(data[,j]))
  }
  return(epmfs)
}


gdtrans <- function(mapped_data){
  # generalized distributional transformation
  transformed_data <- array(dim = dim(mapped_data))
  Fhat_mapped_data <- Fhat(mapped_data)
  pmf_hat_mapped_data <- pmf_hat(mapped_data)
  V <- matrix(runif(n = nrow(mapped_data)*ncol(mapped_data)), ncol = ncol(mapped_data)) 
  for(j in 1:ncol(mapped_data)){
    transformed_data[,j] <- Fhat_mapped_data[[j]](mapped_data[,j]) - pmf_hat_mapped_data[[j]](mapped_data[,j]) + V[,j] * pmf_hat_mapped_data[[j]](mapped_data[,j])
  }
  return(transformed_data)
}

#' G-trasformation function
#' @description function to transform discrete or mixture of discrete and continuous datasets to continuous datasets with marginal normal(0,1). 
#' 
#' @param data The dataset to be transforms. The dataset can be discretein all columns, continuous in all columns or a mixture of continuous columns and discrete columns.
#' @return A transformed continuous dataset with the same copula as the input dataset and margianl normal(0,1).
#' @export
Gtrans <- function(data){
  if (!all(sapply(data, is.numeric))){
    stop("All columns/attributes of the data need to be numeric")
  }
  data_c <- gdtrans(data)
  data_normal <- data.frame(qnorm(data_c, mean = 0, sd = 1))
  colnames(data_normal) <- colnames(data)
  return(data_normal)
}