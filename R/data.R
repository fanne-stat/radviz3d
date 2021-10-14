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
"ceramic"

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

