#' @import rgl
# anchor points coordinates
anchors_sphere <- function(p) {
    phi <- (sqrt(5) + 1)/2
    anchor_coord <- list()
    anchor_coord[[4]] <- matrix(c(1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1), byrow = T, 
        ncol = 3)/sqrt(3)
    anchor_coord[[6]] <- matrix(c(1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, 
        -1), byrow = T, ncol = 3)
    anchor_coord[[8]] <- matrix(c(1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1, 
        1, -1, 1, -1, -1, -1, -1, -1), byrow = T, ncol = 3)/sqrt(3)
    anchor_coord[[12]] <- matrix(c(0, 1, phi, 0, 1, -phi, 0, -1, phi, 0, -1, -phi, 1, 
        phi, 0, -1, phi, 0, 1, -phi, 0, -1, -phi, 0, phi, 0, 1, -phi, 0, 1, phi, 0, -1, 
        -phi, 0, -1), byrow = T, ncol = 3)/sqrt(1 + phi^2)
    anchor_coord[[20]] <- matrix(c(1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1, 1, 1, -1, 
        -1, 1, -1, 1, -1, -1, -1, -1, -1, 0, 1/phi, phi, 0, 1/phi, -phi, 0, -1/phi, phi, 
        0, -1/phi, -phi, 1/phi, phi, 0, -1/phi, phi, 0, 1/phi, -phi, 0, -1/phi, -phi, 
        0, phi, 0, 1/phi, -phi, 0, 1/phi, phi, 0, -1/phi, -phi, 0, -1/phi), byrow = T, 
        ncol = 3)/sqrt(3)
    
    Fibonacci_set <- function(p) {
        anchor <- matrix(nrow = p, ncol = 3)
        for (i in 1:p) {
            anchor[i, 3] <- (2 * i - 1)/p - 1
            anchor[i, 1] <- sqrt(1 - anchor[i, 3]^2) * cos(2 * pi * i * phi^{
                -1
            })
            anchor[i, 2] <- sqrt(1 - anchor[i, 3]^2) * sin(2 * pi * i * phi^{
                -1
            })
        }
        return(anchor)
    }
    
    if (p %in% c(4, 6, 8, 12, 20)) {
        anchors <- anchor_coord[[p]]
    } else {
        anchors <- Fibonacci_set(p)
    }
    return(anchors)
}


# transform data for radial visualization
radial_tranform <- function(data, anchor) {
    data <- as.matrix(data)
    data <- apply(data, MARGIN = 2, FUN = function(x) (x - min(x))/(max(x) - min(x)))
    p <- ncol(data)
    S <- apply(data, MARGIN = 1, sum)
    nominator <- data %*% anchor
    denominator <- matrix(rep(S, each = 3), byrow = T, ncol = 3)
    return(nominator/denominator)
}


# functions to find the location to put labels for class
farthest.point.from.k.centers <- function(x, ctr) {
    
    x.arr <- array(x, dim = c(dim(x), nrow(ctr)))
    ctr.arr <- array(rep(ctr, each = nrow(x)), dim = dim(x.arr))
    
    each.dist <- apply(X = (x.arr - ctr.arr)^2, MARGIN = c(1, 3), FUN = sum)
    
    ## find closest distance of each point to a center
    
    x.min <- apply(X = each.dist, MARGIN = 1, FUN = min)
    
    ## find point that is the fathest from any center
    
    which.max(x.min)
}

separated.class.points <- function(x, cl) {
    ## find most separated points that belong to the different groups
    
    tmp.arr <- apply(X = x, MARGIN = 1, FUN = function(x) (sum(x^2)))
    
    id <- which.max(tmp.arr)
    
    x.arr <- x
    inc.arr <- matrix(x[id, ], ncol = ncol(x))
    
    cls <- cl[id]
    tmp.cl <- cl
    
    for (i in 1:length(unique(cl))) {
        x.arr <- x.arr[tmp.cl != tmp.cl[id], ]
        tmp.cl <- tmp.cl[tmp.cl != tmp.cl[id]]
        id <- farthest.point.from.k.centers(x = x.arr, ctr = inc.arr)
        cls <- c(cls, tmp.cl[id])
        inc.arr <- rbind(inc.arr, x.arr[id, ])
    }
    list(cls, inc.arr)
}

# functions for optimal order of anchor points
estimate.pars <- function(x, cl) {
    p <- ncol(x)
    id <- as.integer(cl)
    K <- max(id)
    ## estimate mixture parameters for calculation of overlap
    Pi <- prop.table(tabulate(id))
    Mu <- t(sapply(1:K, function(k) {
        colMeans(x[id == k, ])
    }))
    S <- sapply(1:K, function(k) {
        var(x[id == k, ])
    })
    dim(S) <- c(p, p, K)
    list(Pi = Pi, Mu = Mu, S = S)
}

prmtd.pars <- function(idx, parlist) {
    Mu <- parlist$Mu[, idx]
    S <- parlist$S[idx, idx, ]
    list(Mu = Mu, S = S)
}

pars.3d <- function(parlist, projmat) {
    Mu <- parlist$Mu %*% projmat
    S <- array(dim = c(3, 3, dim(parlist$S)[3]))
    for (i in 1:(dim(S)[3])) S[, , i] <- t(projmat) %*% parlist$S[, , i] %*% projmat
    list(Mu = Mu, S = S)
}

frobenius.distance <- function(mat1, mat2) (sqrt(sum(diag((mat1 - mat2) %*% t(mat1 - 
    mat2)))))

optimal_3d_anchor_order <- function(x, cl, proj.mat) {
    ll <- estimate.pars(x, cl)
    target.overlap <- MixSim::overlap(Pi = ll$Pi, Mu = ll$Mu, S = ll$S)
    
    x.minmax <- apply(x, MARGIN = 2, FUN = function(x) (x - min(x))/(max(x) - min(x)))
    y <- t(apply(X = x.minmax, MARGIN = 1, FUN = function(z) (z/sum(z))))
    lly <- estimate.pars(y, cl)
    
    p <- ncol(x)
    # (p-1)! permutations. Keep the last index fixed due to the symmetry of sphere
    list.idx <- cbind(gtools::permutations(p - 1, p - 1), p)
    curr.idx <- NULL
    curr.mse <- Inf
    
    for (i in 1:nrow(list.idx)) {
        idx <- list.idx[i, ]
        ll.idx <- prmtd.pars(idx = idx, parlist = lly)
        ll.3d <- pars.3d(parlist = ll.idx, projmat = proj.mat)
        idx.overlap <- MixSim::overlap(Pi = ll$Pi, Mu = ll.3d$Mu, S = ll.3d$S)
        
        idx.mse <- frobenius.distance(idx.overlap$OmegaMap, target.overlap$OmegaMap)
        if (idx.mse < curr.mse) {
            curr.idx <- idx
            curr.mse <- idx.mse
        }
    }
    curr.idx
}


sqrt_scale <- function(x){
  x/sqrt(sum(x^2)) * sqrt(sqrt(sum(x^2)))
}


#' 3D Radial Visualization function
#'
#' @param data The dataset to visualize. Each row is an observation.
#' @param cl The class identification for each observation. The length of \code{cl} should be the same as the number of rows of \code{data}. If specified, different classes would be visualized with different colors.
#' @param domrp Logical. If true, MRP is applied to the origianl dataset. The default number of PCs used is \code{npc = 4}.
#' @param doGrans Logical. If true, Gtrans is applied to the origianl dataset. @seealso \code{\link{Gtrans}}.
#' @param sqrt_scale Logical. If true, the distance of the points to be visualization will be augmented to squre root of the orginal distance to make points further away from the origin.
#' @param color The colors for different classes. If not specified, \code{rainbow} is used.
#' @param colorblind Logical.The colors for different classes.If true, poits are colorblind friendly.If false, \code{rainbow} is used.
#' @param axis Logical.If true, Cartesian axes would be plotted.
#' @param pradius The radius of the data point in RadViz3D. The default value is 0.01.
#' @param with.coord.labels Logical. If true, labels of coordinates will be added to the visualization.
#' @param coord.labels The labels for components of the dataset. When \code{domrp = TRUE}, the coord.labels will be changed to "Xi" representing the the ith direction obtained with MRP.
#' @param coord.font The font for labels of components.
#' @param coord.cex The size of the labels of components.
#' @param with.class.labels Logical. If true, class labels will be added to the visualization.
#' @param class.labels The labels for different classes in the dataset.
#' @param class.labels.locations Locations to put labels for each class. If not specified, an optimal location for each class would be calculated.
#' @param opt.anchor.order Logical. If true, the optimal order of anchor points corresponding to the components would be calculated. This is a very time consuming procedure. Not recommended if the number of components is larger then 6.
#' @param ... Some other parameters from \link{mrp} and \link{Gtrans}.
#' @return A list with the elements
#' \item{mrp.res}{The result of MRP is the argument \code{domrp = TRUE}. See also \code{\link{mrp}}.}
#' @examples
#' radialvis3d(data = iris[,-5], cl = iris[,5], domrp = T)
#' @export
radialvis3d <- function(data, domrp = T, doGtrans = F, sqrt_scale=F, cl = NULL, color = NULL, colorblind = FALSE, axis = FALSE, pradius = 0.01, with.coord.labels = T, coord.labels = NULL, coord.font = 2, coord.cex = 1.1, with.class.labels = T,
    class.labels = levels(factor(cl)), class.labels.locations = NULL, opt.anchor.order = FALSE, ...) {
    if (is.null(cl)) {
        cl <- as.factor(1)
    } else {
        cl <- as.factor(cl)
    }
   
    class <- levels(cl)
    
    if (is.null(color)) {
        if (colorblind == TRUE) {
            color <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7")
        }
        if (colorblind == FALSE) {
            color <- rainbow(length(class))
        }
    }
    
    if (length(color) < length(class)) {
        stop("Not enough colors supplied\n")
    }
    
    res <- list()
    
    # browser()
    if (doGtrans){
      data <- Gtrans(data, cl = cl, ...)
    }

    if (domrp){
      if (length(levels(cl)) <= 1){
        stop("Cannot do MRP with no class information or only one class\n")
      }
      mrp_res <- mrp(data = data, cl = cl, ...)
      res$mrp.res <- mrp_res
      data <- mrp_res$projected_df
      coord.labels <- NULL
      for (i in colnames(data)){
        coord.labels <- c(coord.labels, parse(text = paste0("italic(", i, ")")))
      }
    }
    
    anchors <- anchors_sphere(ncol(data))
    
    # optimal order of anchor points
    idx.opt <- 1:ncol(data)
    if ((length(class) > 1) & opt.anchor.order) {
        idx.opt <- optimal_3d_anchor_order(x = data, cl = cl, proj.mat = anchors)
    }
    
    data_trans <- radial_tranform(data[, idx.opt], anchors)
    
    # sqrrt of each coordiantes
    if (sqrt_scale) data_trans <- t(apply(data_trans, 1, sqrt_scale))
    
    max_distance <- max(apply(data_trans, MARGIN = 1, FUN = function(x) sqrt(sum(x^2))))
    
    radius <- 1.05 * max_distance
    
    # Create rgl plot
    rgl.open()
    rgl.bg(color = "white", alpha = 0.01)
    for (i in 1:length(class)) {
        rgl.spheres(data_trans[cl == class[i], ], r = pradius, color = color[i])
    }
    for (p in 1:ncol(data)) {
        rgl.lines(rbind(rep(0, 3), radius * anchors[p, ]), col = "gray40")
    }
    rgl.spheres(x = 0, y = 0, z = 0, r = radius, color = "grey", add = TRUE, alpha = 0.02)
    rgl.points(radius * anchors, color = "black")

    
    # add labels for coordinates
    if (with.coord.labels){
      if (is.null(coord.labels)){
        coord.labels <- colnames(data)
      }
      text3d(1.1 * radius * anchors, texts = coord.labels, font = coord.font , cex = coord.cex)
    }
    

    # If axis should be plotted
    if (axis == TRUE) {
        # Add axes
        rgl.lines(c(-1, 1), c(0, 0), c(0, 0), color = "black")
        rgl.lines(c(0, 0), c(-1, 1), c(0, 0), color = "black")
        rgl.lines(c(0, 0), c(0, 0), c(-1, 1), color = "black")
    }
    
    # add labels for classes

    if (length(levels(cl)) > 1) {
      if (with.class.labels){
        if (is.null(class.labels.locations)) {
            ll <- separated.class.points(x = data_trans, cl = as.numeric(factor(cl)))

            text3d(ll[[2]] + c(-0.01, -0.01, -0.01), color = color[ll[[1]]], texts = class.labels[ll[[1]]], ...)
        } else text3d(class.labels.locations, color = color[ll[[1]]], texts = class.labels[ll[[1]]], ...)
    }}

    view3d(theta = 60, phi = 60)
    # return(res)
}

