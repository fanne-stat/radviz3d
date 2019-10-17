Fhat <- function(data){
  # returns a list of emperical cdfs. The list is of the same length of ncol(data)
  ecdfs <- list()
  
  for (j in 1:ncol(data)){
    ecdfs[[j]] <- ecdf(data[,j])
  }
  return(ecdfs)
}

gdtrans <- function(mapped_data){
  # generalized distributional transformation
  transformed_data <- array(dim = dim(mapped_data))
  n <- dim(mapped_data)[1]
  Fhat_mapped_data <- Fhat(mapped_data)
  V <- matrix(runif(n = nrow(mapped_data)*ncol(mapped_data)), ncol = ncol(mapped_data)) 
  for(j in 1:ncol(mapped_data)){
    sdata <- sort(mapped_data[,j])
    pmf <- tabulate(match(sdata, unique(sdata))) / n
    sindex <- match(mapped_data[,j], unique(sdata))
    transformed_data[,j] <- Fhat_mapped_data[[j]](mapped_data[,j]) - pmf[sindex] + V[,j] * pmf[sindex]
  }
  return(transformed_data)
}

#' G-trasformation function
#' @description function to transform discrete or mixture of discrete and continuous datasets to continuous datasets with marginal normal(0,1). 
#' 
#' @param data The dataset to be transforms. The dataset can be discretein all columns, continuous in all columns or a mixture of continuous columns and discrete columns.
#' @param cl The class information of the dataset. This is not required when \code{VariableSelection = FALSE}.
#' @param VariableSelection Logical. If true, anova will be performed to each variable to see whether there is a difference among groups for that variable. The varaible associated with Bonferroni adjusted p-value larger than a threshold will be removed.
#' @param p_threshold The threshold for adjusted p-value in variable selection when \code{VariableSelection = TRUE}.
#' @return A transformed continuous dataset with the same copula as the input dataset and margianl normal(0,1).
#' @export
Gtrans <- function(data, cl = NULL ,VariableSelection = FALSE, p_threshold = 0.05, ...){
  if (!all(sapply(data, is.numeric))){
    stop("All columns/attributes of the data need to be numeric")
  }
  data_c <- gdtrans(data)
  data_normal <- data.frame(qnorm(data_c, mean = 0, sd = 1))
  colnames(data_normal) <- colnames(data)
  
  if(VariableSelection){
    if(is.null(class)){
      stop("class information required for variable selection!")
    }
    mat.lm.pval <- apply(data_normal, 2, function(x, gr)(anova(lm(x~gr))[5][[1]][1]), gr = as.factor(cl))
    data_normal <- data_normal[,p.adjust(mat.lm.pval, method = "BH") < p_threshold]
  }
  return(data_normal)
}
