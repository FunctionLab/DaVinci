# https://github.com/jmonlong/Hippocamplus/blob/master/content/post/2018-06-09-ClusterEqualSize.Rmd


#' Tile the slice for scalable implementation
#' @export
tile_the_slice <- function(coord, random.seed = 1, L2_number = 1000){
  if (is.null(rownames(coord))){
    stop("Input doesn't have rownames.")
  }#if

  #pre-set
  L1_size <- 7000
  #L2_number <- 1000 #about 1000 spots in the problem
  message("L2_number is not recommended to be higher than 2000. Higher the number is, the more granularity you expect.")

  #derived
  L1_number <- nrow(coord)/L1_size #number
  L2_size <- nrow(coord) %/% L2_number #balanced.cluster.size
  message(paste0("Tile size is ", L2_size))


  #calculate the L1 grid parameter
  #################################
  bbox <- sf::st_bbox(c(xmin = min(coord[,1]), ymin = min(coord[,2]), xmax = max(coord[,1]), ymax = max(coord[,2]) ))
  xy_ratio <- (bbox$xmax-bbox$xmin)/(bbox$ymax-bbox$ymin)
  ny <- ceiling(sqrt(L1_number/xy_ratio))
  nx <- ceiling(ny*xy_ratio)
  grid <- sf::st_make_grid(bbox, n = c(nx, ny), square = T)
  pos <- sf::st_as_sf(as.data.frame(coord), coords = c("array_row","array_col"))

  ids_L1 <- unlist(lapply(sf::st_intersects(pos, grid), function(sublist) sublist[1]))
  names(ids_L1) <- rownames(coord)


  #L2 fine grid
  #################################
  ids_L2 <- ids_L1
  ids_L1_unique <- sort(unique(ids_L1))

  for (ii in 1:length(ids_L1_unique)){
    temp <- which(ids_L1==ids_L1_unique[ii])
    if (length(temp)<L2_size){
      ids_L2[temp] <- 1
    }else{
      ids_L2[temp] <- swk(coord[temp,], method = "balanced", L2norm = F, balanced.cluster.size = L2_size, random.seed = random.seed)
    }#else

  }#for ii

  #combine together
  #################################
  grid_ids <- paste0(ids_L1, "_", ids_L2)
  mapping <- 1:length(unique(grid_ids))
  names(mapping) <- unique(grid_ids)
  grid_ids <- mapping[grid_ids]

  names(grid_ids) <- rownames(coord)
  return(grid_ids)
}#tile_the_slice






#' Same Size Clustering
#'
#' This is a wrapper for several implementation that classify samples into
#' same size clusters, the details please see [this blog](http://jmonlong.github.io/Hippocamplus/2018/06/09/cluster-same-size/).
#' The source code is modified based on code from the blog.
#'
#' @param mat a data/distance matrix.
#' @param diss if `TRUE`, treat `mat` as a distance matrix.
#' @param clsize integer, number of sample within a cluster.
#' @param algo algorithm.
#' @param method method.
#'
#' @return a vector.
#' @export
#'
#' @examples
#' set.seed(1234L)
#' x <- rbind(
#'   matrix(rnorm(100, sd = 0.3), ncol = 2),
#'   matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2)
#' )
#' colnames(x) <- c("x", "y")
#'
#' y1 <- same_size_clustering(x, clsize = 10)
#' y11 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, diss = TRUE)
#'
#' y2 <- same_size_clustering(x, clsize = 10, algo = "hcbottom", method = "ward.D")
#'
#' y3 <- same_size_clustering(x, clsize = 10, algo = "kmvar")
#' y33 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, algo = "kmvar", diss = TRUE)
same_size_clustering <- function(mat, diss = FALSE, clsize = NULL,
                                 algo = c("nnit", "hcbottom", "kmvar"),
                                 method = c(
                                   "maxd", "random", "mind", "elki",
                                   "ward.D", "average", "complete", "single"
                                 )) {
  stopifnot(is.numeric(clsize))

  algo <- match.arg(algo)
  method <- match.arg(method)
  do.call(algo, args = list(mat = mat, diss = diss, clsize = clsize, method = method))
}

nnit <- function(mat,
                 clsize = NULL,
                 diss = FALSE,
                 method = "maxd") {
  stopifnot(is.logical(diss))

  clsize.rle <- rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat) / clsize))))
  clsize <- clsize.rle$lengths
  lab <- rep(NA, nrow(mat))
  if (isFALSE(diss)) {
    dmat <- as.matrix(dist(mat))
  } else {
    dmat <- mat
  }
  cpt <- 1
  while (sum(is.na(lab)) > 0) {
    lab.ii <- which(is.na(lab))
    dmat.m <- dmat[lab.ii, lab.ii]
    ii <- switch(method,
                 maxd = which.max(rowSums(dmat.m)),
                 mind = which.min(rowSums(dmat.m)),
                 random = sample.int(nrow(dmat.m), 1),
                 stop("unsupported method in 'nnit'!")
    )
    lab.m <- rep(NA, length(lab.ii))
    lab.m[head(order(dmat.m[ii, ]), clsize[cpt])] <- cpt
    lab[lab.ii] <- lab.m
    cpt <- cpt + 1
  }
  if (any(is.na(lab))) {
    lab[which(is.na(lab))] <- cpt
  }
  lab
}#same_size_clustering


kmvar <- function(mat,
                  clsize = NULL,
                  diss = FALSE,
                  method = "maxd") {
  stopifnot(is.logical(diss))

  k <- ceiling(nrow(mat) / clsize)
  if (isFALSE(diss)) {
    km.o <- kmeans(mat, k)
    # distance to centers
    centd <- lapply(1:k, function(kk) {
      euc <- t(mat) - km.o$centers[kk, ]
      sqrt(apply(euc, 2, function(x) sum(x^2)))
    })
    centd <- matrix(unlist(centd), ncol = k)
  } else {
    message("PAM algorithm is applied when input distance matrix.")
    pam.o <- cluster::pam(mat, k, diss = TRUE)
    # medoids
    # distance to medoids
    centd <- mat[, pam.o$id.med, drop = FALSE]
  }

  labs <- rep(NA, nrow(mat))
  clsizes <- rep(0, k)

  ptord <- switch(method,
                  maxd = order(-apply(centd, 1, max)),
                  mind = order(apply(centd, 1, min)),
                  random = sample.int(nrow(mat)),
                  elki = order(apply(centd, 1, min) - apply(centd, 1, max)),
                  stop("unsupported method in 'kmvar'!")
  )

  for (ii in ptord) {
    bestcl <- which.max(centd[ii, ])
    labs[ii] <- bestcl
    clsizes[bestcl] <- clsizes[bestcl] + 1
    if (clsizes[bestcl] >= clsize) {
      centd[, bestcl] <- NA
    }
  }
  return(labs)
}#kmvar


#' @importFrom stats as.dist
hcbottom <- function(mat,
                     clsize = NULL,
                     diss = FALSE,
                     method = "ward.D") {
  stopifnot(is.logical(diss))

  method <- match.arg(method, choices = c("ward.D", "average", "complete", "single"))
  if (isFALSE(diss)) {
    dmat <- as.matrix(dist(mat))
  } else {
    dmat <- mat
  }
  clsize.rle <- rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat) / clsize))))
  clsizes <- clsize.rle$lengths
  cpt <- 1
  lab <- rep(NA, nrow(mat))
  for (clss in clsizes[-1]) {
    lab.ii <- which(is.na(lab))
    hc.o <- hclust(as.dist(dmat[lab.ii, lab.ii]), method = method)
    clt <- 0
    ct <- length(lab.ii) - clss
    while (max(clt) < clss) {
      cls <- cutree(hc.o, ct)
      clt <- table(cls)
      ct <- ct - 1
    }
    cl.sel <- which(cls == as.numeric(names(clt)[which.max(clt)]))
    lab[lab.ii[head(cl.sel, clss)]] <- cpt
    cpt <- cpt + 1
  }
  lab[is.na(lab)] <- cpt
  lab
}#hcbottom





