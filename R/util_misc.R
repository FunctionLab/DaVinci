require(mclust)

#source("/mnt/ceph/users/wmao1/KPMP/code/PLIER_Ahmed/speed_sylvester/L4_adaptive/scaling_1/utils_misc.R")



CM <- function(x,y){
  include.index <- which(!is.na(x) & !is.na(y))
  return(aricode::ARI(x[include.index], y[include.index]))
}#CM



L2Norm <- function(mat, MARGIN = 1){
  normalized <- sweep(
    x = mat,
    MARGIN = MARGIN,
    STATS = apply(
      X = mat,
      MARGIN = MARGIN,
      FUN = function(x){
        sqrt(x = sum(x ^ 2))
      }
    ),
    FUN = "/"
  )
  normalized[!is.finite(x = normalized)] <- 0
  return(normalized)
}#L2Norm




#wrapper function with different clustering algorithm as support
#partition <- swk(ICAp.res$B, method = "mclust", L2norm = T, mclust.model = "EEE", mclust.num = gt.number[pid], random.seed = 1)
#partition <- swk(ICAp.res$B, method = "kmeans", L2norm = T, km.num = gt.number[pid], random.seed = 1)
#partition <- swk(ICAp.res$B, method = "leiden", L2norm = T,  ld.resolution = 0.2, random.seed = 1)
#https://github.com/scverse/scanpy/issues/1531

swk <- function(kmeans.input, method = "mclust", L2norm = T, zscore = F, mclust.model = "EEE", mclust.num = NULL, km.num = NULL, ld.num.neighbors = 20, ld.resolution = NULL, ld.NN = "SNN", balanced.num = NULL, balanced.cluster.size = NULL, random.seed = 1, weights = NULL){

  #flip the coordinates
  if (ncol(kmeans.input) > nrow(kmeans.input) ){
    kmeans.input <- t(kmeans.input)
  }#if

  #normalization
  #####################################
  if (L2norm){
    kmeans.input <- L2Norm(kmeans.input, MARGIN = 2)
    #kmeans.input <- L2Norm(kmeans.input, MARGIN = 1)
  }else if (zscore){
    kmeans.input <- scale(kmeans.input)
  }#else

  #addd weights
  if (!is.null(weights)){
    if (length(weights)!=ncol(kmeans.input)){
      stop("weights don't match the dimension")
    }else{
      kmeans.input <- sweep(kmeans.input, 2, weights, FUN = "*")

    }#else
  }#if


  #clustering step
  #####################################
  set.seed(random.seed)

  if (method == "mclust"){

    if (is.null(mclust.num)){
      stop("Number of clusters is not set in mclust")
    }#if

    mclust.res <- mclust::Mclust(kmeans.input, G = mclust.num, modelNames = mclust.model, verbose = F)
    partition <- mclust.res$classification


    #if is.null
    if (is.null(partition)){
      count <- 0
      message("Repeat mclust")
      while ( (is.null(partition))& (count < 3)) {
        mclust.res <- mclust::Mclust(kmeans.input, G = mclust.num, modelNames = mclust.model)
        partition <- mclust.res$classification
        count <- count+1
      }#while
    }#if

    #if still null
    ##############
    if (is.null(partition)){

      message("Supervised mclust")
      set.seed(1)

      training.index <- sample(1:nrow(kmeans.input), nrow(kmeans.input)*0.2)
      test.index <- setdiff(1:nrow(kmeans.input), training.index)
      input.training <- kmeans.input[training.index,]
      input.test <- kmeans.input[test.index,]

      mclust.res <- mclust::Mclust(input.training, G = mclust.num, modelNames = mclust.model, verbose = F)
      p.training <- mclust.res$classification

      test.res <- predict(mclust.res, newdata = input.test)
      p.test <- test.res$classification

      partition <- rep(NA, nrow(kmeans.input))
      partition[training.index] <- p.training
      partition[test.index] <- p.test
    }#if


  }else if (method == "kmeans"){

    if (is.null(km.num)){
      stop("Number of clusters is not set in kmeans")
    }#if

    km.res <- kmeans(kmeans.input, centers = km.num, nstart = 25)
    partition <- km.res$cluster

  }else if (method == "leiden"){

    message("Number of neighbors is ", ld.num.neighbors)

    if (is.null(ld.resolution)){
      stop("Resolution parameter is not set in leiden")
    }#if


    snn.res <- Seurat::FindNeighbors(kmeans.input, k.param = ld.num.neighbors, return.neighbor = F, compute.SNN = T, verbose = F)
    if (ld.NN =="SNN"){
      message("SNN")
      leiden.input <- snn.res$snn

    }else if (ld.NN=="SNN.binary"){
      message("SNN.binary")
      leiden.input <- snn.res$snn
      leiden.input[leiden.input!=0] <- 1

    }else if (ld.NN=="knn"){
      message("kNN")
      leiden.input.index <- snn.res$nn
      leiden.input <- as.matrix(dist(kmeans.input))

      #apply gaussian kernel to convert to similarity score
      sigma <- 1
      leiden.input <- exp(-leiden.input^2/(2*sigma^2))
      leiden.input[which(as.matrix(leiden.input.index)==0)] <- 0

    }else if (ld.NN=="knn.binary"){
      message("kNN.binary")
      leiden.input <- snn.res$nn

    }#else if

    partition <- leiden::leiden(leiden.input, resolution_parameter = ld.resolution)

  }else if(method == "balanced"){

    if (is.null(balanced.cluster.size)){
      if ( is.null(balanced.num) ){
        stop("Cannot determine the cluster size")
      }else{
        balanced.cluster.size <- nrow(kmeans.input) %/% balanced.num
      }#else
    }#if

    balanced.res <- same_size_clustering(kmeans.input, diss = F, clsize = balanced.cluster.size, algo = "nnit", method = "maxd")
    partition <- balanced.res

  }else{
    stop("The method is not included in the pre-defined list")
  }#else

  #return the clusterings
  #####################################
  if (is.null(partition)){
    stop("Number of unique clusters is 0")
  }else{
    print(paste0("Number of unique clusters is ", length(unique(partition))))

    names(partition) <- rownames(kmeans.input)
    return(partition)
  }#else

}#swk






normF <- function(x){
  return(sqrt(sum(x^2)))
}#normF



soft_threshold <- function(x, t){
  return(sign(x)*sapply(abs(x)-t, function(x){max(x,0)}))
}#soft_threshold




flr.batch <- function(data, D=NULL, L, lambda = 1e-2, rho = 10, tol = 1e-2, max.iter = 100, kmeans =F, verbose = F){
  if (ncol(data) < nrow(data)){
    data <- t(data)
  }#if

  for (ii in 1:nrow(data)){
    data[ii, ] <- flr(data[ii,], D=D, L=L, lambda = lambda, rho = rho, tol = tol, max.iter = max.iter, kmeans = kmeans, verbose = verbose)
  }#for ii

  return(data)
}#flr.batch




#https://www.stat.cmu.edu/~ryantibs/convexopt-F15/lectures/21-dual-meth.pdf
#https://web.stanford.edu/~boyd/papers/pdf/network_lasso.pdf

flr <- function(y, D = NULL, L, lambda = 1e-2, rho = 10, tol = 1e-2, max.iter = 100, kmeans = F, verbose = F){

  #construct D based on the Laplacian matrix
  ###############################################
  if (is.null(D)){
    adj <- -L
    diag(adj) <- 0

    #all(rownames(L)==colnames(L))
    #all(rownames(L)==colnames(res1$B))
    grp <- igraph::graph_from_adjacency_matrix(adj)
    edge.list <- igraph::get.edgelist(grp)
    D <- matrix(0, nrow = nrow(edge.list), ncol = ncol(L))
    colnames(D) <- colnames(L)
    for (ii in 1:nrow(edge.list)){
      D[ii, edge.list[ii,1]] <- -1
      D[ii, edge.list[ii,2]] <- 1
    }#for ii
  }#if is.null(D)



  #parameter initalization
  ###############################################
  iter <- 0
  diff <- 1

  y <- matrix(y, ncol = 1)
  beta <- matrix(rep(0, nrow(y)), ncol = 1)

  if (kmeans){
    km.res <- kmeans(y, centers = 2 ,nstart = 25)
    beta[which(km.res$cluster==1),1] <- km.res$centers[1,]
    beta[which(km.res$cluster==2),1] <- km.res$centers[2,]
  }#kmeans

  alpha <- D%*%beta
  w <- matrix(rep(0, nrow(D)), ncol = 1)


  #pre-calculate
  ##################
  Iden <- diag(nrow(y))
  DT <- t(D)
  D.sparse <- as(D, "sparseMatrix")
  DT.sparse <- as(DT, "sparseMatrix")
  #ptm <- proc.time()
  DTD.sparse <- DT.sparse %*% D.sparse #0.006s vs 240s
  #print(proc.time()-ptm)

  #ptm <- proc.time()
  left_inv <- as.matrix(Matrix::solve(rho*DTD.sparse+Iden)) #0.485s vs 35s
  #left_inv <- solve(Iden+rho*DTD)
  #print(proc.time()-ptm)




  #start
  ##################
  while( (diff >= tol)  & (iter < max.iter)){

    if (verbose){
      print(diff)
    }#if

    beta.old <- beta

    #update steps
    #####################
    beta <- left_inv %*% (y+rho*DT%*%(alpha-w))

    Dbeta_w <- D%*%beta+w

    alpha <- soft_threshold(Dbeta_w , lambda/rho)

    w <- Dbeta_w-alpha


    #wrap up
    #####################
    iter <- iter+1
    diff <- normF(beta.old - beta)

  }#while

  return(beta)
}#flr








moranI <- function(x, L){
  #recover the adj matrix from L
  A <- -L
  diag(A) <- 0

  #calculate the moran'I
  #https://rspatial.org/rosu/Chapter7.html

  terra::autocor(x, A, method = "moran")
}#moranI


moranI_all <- function(x, L){

  res <- c()
  if (nrow(x)>ncol(x)){
    x <- t(x)
  }#if
  for (i in 1:nrow(x)){
    res <- c(res, moranI(x[i,], L))
  }#for i
  return(res)
}#moranI_all







#use swk
split_chunk <- function(x, cluster.size, random.seed = 1, reference_set = NULL, reference_set_opt){


  num.of.chunk <- length(x)%/% cluster.size

  set.seed(random.seed)

  if (is.null(reference_set)){
    message("random")
    partition <- sample(1:length(x)%%num.of.chunk)

  }else{

    if (length(x)!=nrow(reference_set)){
       stop("Dimensions don't match")
    }#if


    if (reference_set_opt$method == "leiden"){

      message("leiden")
      partition <- swk(reference_set, method = "leiden", L2norm = reference_set_opt$L2norm, ld.resolution = reference_set_opt$ld.resolution, ld.num.neighbors = reference_set_opt$ld.num.neighbors, ld.NN = reference_set_opt$ld.NN, random.seed = random.seed)

    }else if (reference_set_opt$method == "mclust"){

      message("mclust")
      partition <- swk(reference_set, method = "mclust", L2norm = reference_set_opt$L2norm, mclust.model = reference_set_opt$mclust.model,  mclust.num = num.of.chunk, random.seed = random.seed)

    }else if (reference_set_opt$method == "balanced"){

      message("balanced")
      partition <- swk(reference_set, method = "balanced", L2norm = reference_set_opt$L2norm, balanced.cluster.size = cluster.size, random.seed = random.seed)

    }#else if


  }#else

  to_report <- split(x, partition)
  return(to_report)

}#split_chunk








mean.diff <- function(x, index1, index2){
  return( wilcox.test(x[index1], x[index2])$p.value )
}#mean.diff






#single slice: niche can be used as batch
#multiple slices: data can be concatenated matrix; batch can be slice x niche
#reference_set: spot-by-LV/xandy for clustering purpose
#leiden, mclust, balanced set clustering

#' Differential analysis
#' @export
pseudo.default <- function(data, one_hot_encode, batch = NULL, num.of.pseudo, verbose = T, reference_set = NULL, reference_set_opt = list(method="mclust", L2norm = T, ld.resolution=1, ld.num.neighbors=20, ld.NN="SNN", mclust.model = "EEE")){

  if (ncol(data)!=length(one_hot_encode)){
    stop("Dimension doesn't match")
  }#if

  one_index <- which(one_hot_encode=="1")
  zero_index <- which(one_hot_encode!="1")

  if (is.null(num.of.pseudo)){
    dat.one <- data[,one_index]
    dat.zero <- data[,zero_index]

  }else{

    if (is.null(batch)){

      #get the index for one
      one_list <- split_chunk(one_index, cluster.size = num.of.pseudo, reference_set = reference_set[one_index,], reference_set_opt = reference_set_opt)
      zero_list <- split_chunk(zero_index, cluster.size = num.of.pseudo, reference_set = reference_set[zero_index,], reference_set_opt = reference_set_opt)


    }else{

      batch_one <- unique(batch[one_index])
      batch_zero <- unique(batch[zero_index])

      print(batch_one)
      print(batch_zero)

      #get the index for one
      ################################
      one_list <- list()
      for (i in 1:length(batch_one)){
        temp <- which(one_hot_encode==1 & batch==batch_one[i])
        one_list <- c(one_list, split_chunk(temp, cluster.size = num.of.pseudo, reference_set = reference_set[temp,], reference_set_opt = reference_set_opt) )
      }#for i

      #get the index for zero
      ################################
      zero_list <- list()
      for (i in 1:length(batch_zero)){
        temp <- which(one_hot_encode==0 & batch==batch_zero[i])
        zero_list <- c(zero_list, split_chunk(temp, cluster.size = num.of.pseudo, reference_set = reference_set[temp,], reference_set_opt = reference_set_opt) )
      }#for i

    }#else

    #get the pseudo matrix
    ############################################################
    dat.one <- do.call(cbind, lapply(one_list, function(x, data){ apply(data[,x, drop=F], 1,mean)}, data = data) )
    dat.zero <- do.call(cbind, lapply(zero_list, function(x, data){ apply(data[,x, drop= F], 1,mean)}, data = data) )

  }#else, !is.null(num.of.pseudo)

  if (verbose){
    cat("one: ", ncol(dat.one), ", zero: ", ncol(dat.zero),"\n")
  }#if



  #row-wise wilcoxon test
  ############################################################
  dd <- cbind(dat.one, dat.zero)
  index1 <- 1:ncol(dat.one)
  index2 <- (ncol(dat.one)+1):ncol(dd)
  p.val <- apply(dd,1,mean.diff, index1 = index1, index2 = index2)

  return(p.val)
}#pseudo.default



#' Pathway enrichment analysis via EnrichR
#' @export
enrich_visual <- function(gene.list, gene.all, dbs = c("HuBMAP_ASCTplusB_augmented_2022"), number.of.top=10, pathway = "HuBMAP_ASCTplusB_augmented_2022"){
  
  enrich.res <- enrichr(gene.list, databases = dbs, background = gene.all)
  
  dat.to_plot <- enrich.res[[pathway]]
  dat.to_plot <- dat.to_plot[order(dat.to_plot$Adjusted.P.value,decreasing = F)[1:number.of.top],]
  dat.to_plot$FDR <- -log10(dat.to_plot$Adjusted.P.value)
  dat.to_plot$Term <- stringr::str_wrap(dat.to_plot$Term, width = 50)
  
  p <- ggdotchart(dat.to_plot, x = "Term", y = "FDR",
                  color = "#D55E00",
                  palette = c(  "#DDDDDD","#D55E00"),
                  sorting = "descending",
                  add = "segments",
                  add.params = list(color = "#DDDDDD", size = 2),
                  dot.size = 10,
                  label = NULL,
                  font.label = list(color = "white", size = 9,vjust = 0.5),
                  ggtheme = theme_pubr(base_size = 15),
                  rotate=T
  )+geom_hline(yintercept = -log10(0.05), linetype = 2, color = "black")+xlab("")+ylab("-log10(FDR)")+theme(legend.position = "none", axis.text.x = element_text(size = 20), axis.title.y = element_text(face = "bold", size= 30))
  
  p
}#enrich_visual
