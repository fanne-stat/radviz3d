% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radviz3d.R
\name{radialvis3d}
\alias{radialvis3d}
\title{3D Radial Visualization function}
\usage{
radialvis3d(
  data,
  domrp = T,
  doGtrans = F,
  sqrt_scale = F,
  cl = NULL,
  color = NULL,
  pch = 16,
  colorblind = FALSE,
  axes = FALSE,
  point.cex = 1,
  with.coord.labels = T,
  coord.labels = NULL,
  coord.font = 2,
  coord.cex = 1.1,
  with.class.labels = T,
  class.labels = levels(factor(cl)),
  class.labels.locations = NULL,
  opt.anchor.order = FALSE,
  alpha = 0.02,
  lwd = 1,
  axes.col = "black",
  ret.trans = FALSE,
  ...
)
}
\arguments{
\item{data}{The dataset to visualize. Each row is an observation.}

\item{domrp}{Logical. If true, MRP is applied to the origianl dataset. The default number of PCs used is \code{npc = 4}.}

\item{doGtrans}{Logical. If true, Gtrans is applied to the origianl dataset. @seealso \code{\link{Gtrans}}.}

\item{sqrt_scale}{Logical. If true, the distance of the points to be visualization will be augmented to squre root of the orginal distance to make points further away from the origin.}

\item{cl}{The class identification for each observation. The length of \code{cl} should be the same as the number of rows of \code{data}. If specified, different classes would be visualized with different colors.}

\item{color}{The colors for different classes. If not specified, \code{rainbow} is used.}

\item{pch}{The point character to be used. It is an integer of a vector of integers of the same length of the nrow of the dataset. See \code{\link{points}} for a complete list of characters.}

\item{colorblind}{Logical.The colors for different classes.If true, poits are colorblind friendly.If false, \code{rainbow} is used.}

\item{axes}{Logical.If true, Cartesian axes would be plotted.}

\item{point.cex}{The size of the data point in RadViz3D. The default value is 1.}

\item{with.coord.labels}{Logical. If true, labels of coordinates will be added to the visualization.}

\item{coord.labels}{The labels for components of the dataset. When \code{domrp = TRUE}, the coord.labels will be changed to "Xi" representing the the ith direction obtained with MRP.}

\item{coord.font}{The font for labels of components.}

\item{coord.cex}{The size of the labels of components.}

\item{with.class.labels}{Logical. If true, class labels will be added to the visualization.}

\item{class.labels}{The labels for different classes in the dataset.}

\item{class.labels.locations}{Locations to put labels for each class. If not specified, an optimal location for each class would be calculated.}

\item{opt.anchor.order}{Logical. If true, the optimal order of anchor points corresponding to the components would be calculated. This is a very time consuming procedure. Not recommended if the number of components is larger then 6.}

\item{alpha}{The alpha value that controls the transparency of the sphere in 3d visulization}

\item{lwd}{The line width in the visualization}

\item{axes.col}{Colors of the axes, if needed to be displayed}

\item{ret.trans}{Logical parameter, returns the Radviz3D transformation if TRUE}

\item{...}{Some other parameters from \link{mrp} and \link{Gtrans} and rgl functions.}
}
\value{
A list with the elements
\item{mrp.res}{The result of MRP is the argument \code{domrp = TRUE}. See also \code{\link{mrp}}.}
}
\description{
3D Radial Visualization function
}
\examples{
radialvis3d(data = iris[,-5], cl = iris[,5], domrp = T)
}
