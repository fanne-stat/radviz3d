#' Compositions of ancient Chinese celdon pieces
#' 
#' This dataset contains compositional data of ancient Chinese celdon from Longquan and Jingdezhen kiln from North Song to Ming Dynasties.
#' 
#'@format  A data frame with 19 variables and 88 observations.
#'  \describe{
#'      \item{mf}{Manufacturer of the celdon piece: FLQ for Jingdezhen and LG for Longquan}
#'      \item{era}{The manufacturing time and part of the celdon piece in "time-part" format. There are two different parts (body (b) and glaze (g)) and four times (Song Dynasty (S), Yuan Dynasty (Y), Ming Dynastty(M) and Qing Dynasty (QC)).}
#'      \item{others}{The contents of chemical components.} 
#'  }
#' 
"celadons"

#' Chemical compositions of wine
#' 
#' The dataset contains chemical compositions of wines from 3 cultivars
#' 
#' @format A data frame of 178 observations and 14 variables:
#'     \describe{
#'         \item{cultivar}{The cultivar where the wine is produced}
#'         \item{other variables}{The content of chemical compositions of the wine}
#'     }
"wine"

#' COVID-19 US variants dataset
#' 
#' This is a compositional dataset of the COVID-19 variants in the US from 6/19/2021 to 9/18/2021.
#' 
#' @format A data frame of 140 observations and 14 variables.
#' \describe{
#'     \item{group}{The date.}
#'     \item{type}{weighted}
#'     \item{region}{Region of the US labelled by numbers.}
#'     \item{other variables}{COVID-19 variants compositions.}
#'}
"sarscov2.us.variants"

#' Simulated datasets for testing
#'
#' This is a list containing three simulated datasets, each with 500 observations
#' and 5 classes, used for testing visualization methods.
#'
#' @format A list of 3 data frames, each with 500 observations and 6 variables:
#' \describe{
#'     \item{class}{Factor with 5 levels representing different classes}
#'     \item{X1, X2, X3, X4, X5}{Numeric variables with simulated data}
#' }
"sim_data"

#' Overlap matrices for simulated data
#'
#' This is a list containing three overlap matrices corresponding to the
#' sim_data datasets, showing class separability.
#'
#' @format A list of 3 matrices, each 5x5, representing overlap between classes
"overlap_mat_sim"

