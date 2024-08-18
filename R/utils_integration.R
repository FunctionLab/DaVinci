#reference find the neighbor
#https://www.bioconductor.org/packages/release/bioc/vignettes/BiocNeighbors/inst/doc/approx.html





#' Integrate two modalities
#' @export
BiModalIntegration <- function(modal.1, modal.2, mat.1=modal.1, mat.2=modal.2, n_neighbors = 20, n_neighbors_large = 200, sigma.idx = n_neighbors, snn.far.nn = T, L2norm = T, sd.scale = 1, cross.contant = 1e-4, prune.SNN = 0, kernel.power = 1){
  
  #L2 normalization
  #################################
  if (L2norm){
    if (F){
    
      message("Sample-wise L2Norm")
      modal.1 <- L2Norm(modal.1, MARGIN = 1)
      modal.2 <- L2Norm(modal.2, MARGIN = 1) 
      mat.1 <- L2Norm(mat.1, MARGIN = 1)
      mat.2 <- L2Norm(mat.2, MARGIN = 1)
      
    }else{
      
      message("Column-wise L2Norm")
      modal.1 <- L2Norm(modal.1, MARGIN = 2)
      modal.2 <- L2Norm(modal.2, MARGIN = 2) 
      mat.1 <- L2Norm(mat.1, MARGIN = 2)
      mat.2 <- L2Norm(mat.2, MARGIN = 2)  
      
    }#else
    
  }#if 

  
  #build up all graph
  #################################
  prebuild.index.modal.1 <- BiocNeighbors::buildIndex(modal.1, BNPARAM = BiocNeighbors::AnnoyParam() )
  modal.1.neighbor <- FindNN(modal.1, number.of.NN = n_neighbors, prebuild.index = prebuild.index.modal.1)
  modal.1.neighbor.large <- FindNN(modal.1, number.of.NN = n_neighbors_large, prebuild.index = prebuild.index.modal.1)
  
  
  
  prebuild.index.modal.2 <- BiocNeighbors::buildIndex(modal.2, BNPARAM = BiocNeighbors::AnnoyParam() )
  modal.2.neighbor <- FindNN(modal.2, number.of.NN = n_neighbors, prebuild.index = prebuild.index.modal.2)
  modal.2.neighbor.large <- FindNN(modal.2, number.of.NN = n_neighbors_large, prebuild.index = prebuild.index.modal.2)
  
  
  message("Calculating the modality weights")
  ########################################################################################
  
  #calculate the within and cross-modality distance
  #################################
  
  #nearest distance
  NNdist1 <- modal.1.neighbor$distance[,1]
  NNdist2 <- modal.2.neighbor$distance[,1]
  
  
  #within modality prediction
  mat.1.pred.NN1 <- PredictAssay(mat.1, modal.1.neighbor$index)
  mat.2.pred.NN2 <- PredictAssay(mat.2, modal.2.neighbor$index)
  
  #cross modality prediction
  mat.1.pred.NN2 <- PredictAssay(mat.1, modal.2.neighbor$index)
  mat.2.pred.NN1 <- PredictAssay(mat.2, modal.1.neighbor$index)
  
  
  
  #calculate the distance, within the modality
  impute.dist.1.NN1 <- ImputeDist(mat.1, mat.1.pred.NN1, NNdist1)
  impute.dist.2.NN2 <- ImputeDist(mat.2, mat.2.pred.NN2, NNdist2)
  
  #calculate the distance, across the modality
  impute.dist.1.NN2 <- ImputeDist(mat.1, mat.1.pred.NN2, NNdist1)
  impute.dist.2.NN1 <- ImputeDist(mat.2, mat.2.pred.NN1, NNdist2)
  
  
  
  #calculate the kernel width
  #################################
  if (snn.far.nn){
    snn.matrix.1 <- Seurat:::ComputeSNN(nn_ranked = modal.1.neighbor$index, prune = prune.SNN)
    snn.matrix.2 <- Seurat:::ComputeSNN(nn_ranked = modal.2.neighbor$index, prune = prune.SNN)
    
    snn.width.1 <- Seurat:::ComputeSNNwidth(snn.graph = snn.matrix.1, k.nn = n_neighbors, l2.norm = F, embeddings = mat.1, nearest.dist = NNdist1)
    
    snn.width.2 <- Seurat:::ComputeSNNwidth(snn.graph = snn.matrix.2, k.nn = n_neighbors, l2.norm = F, embeddings = mat.2, nearest.dist = NNdist2)
    
    sd.1 <- snn.width.1*sd.scale
    sd.2 <- snn.width.2*sd.scale
    
    
  }else{
    
    #calculate based on the sigma.idx neighbor
    sd.1 <- (modal.1.neighbor$distance[,sigma.idx]-NNdist1)*sd.scale #wrong implementation has been corrected here
    sd.2 <- (modal.2.neighbor$distance[,sigma.idx]-NNdist2)*sd.scale
    
  }#else
  
  #calculate within and cross modality kernel, and modality weights
  #################################
  
  
  kernel.1.NN1 <- exp(-1*impute.dist.1.NN1/sd.1)
  kernel.2.NN2 <- exp(-1*impute.dist.2.NN2/sd.2)
  
  kernel.1.NN2 <- exp(-1*impute.dist.1.NN2/sd.1)
  kernel.2.NN1 <- exp(-1*impute.dist.2.NN1/sd.2)
  
  
  #compute the modality weights
  score.1 <- kernel.1.NN1/(kernel.1.NN2+cross.contant)
  score.2 <- kernel.2.NN2/(kernel.2.NN1+cross.contant) #typo here has been corrected
  score.1 <- MinMax(score.1, val.min = 0, val.max = 200)
  score.2 <- MinMax(score.2, val.min = 0, val.max = 200)
  
  score.sum <- exp(score.1)+exp(score.2)
  weight.1 <- exp(score.1)/score.sum
  weight.2 <- exp(score.2)/score.sum
  
  
  
  message("Calculating the weighted KNN and SNN")
  ########################################################################################
  
  
  #figure out the union set of neighbors - n_neighbors_large
  modal.union.neighbor.index <- list()
  for (ii in 1:nrow(modal.1.neighbor.large$index)){
    modal.union.neighbor.index[[ii]] <- union(modal.1.neighbor.large$index[ii,], modal.2.neighbor.large$index[ii,])
  }#for ii
  
  
  
  #calculate the new distance matrix in each modality given the union set
  modal.1.union.dist <- neighbor_distance_calculation(modal.union.neighbor.index, modal.1, nearest.dist = NNdist1)
  modal.2.union.dist <- neighbor_distance_calculation(modal.union.neighbor.index, modal.2, nearest.dist = NNdist2)
  
  
  
  #weighted by the modal weights, and then take the total sum
  modal.1.union.dist.weighted <- lapply(1:length(modal.1.union.dist), function(x){ exp(-1*(modal.1.union.dist[[x]]/sd.1[x])** kernel.power)*weight.1[x]})
  modal.2.union.dist.weighted <- lapply(1:length(modal.2.union.dist), function(x){ exp(-1*(modal.2.union.dist[[x]]/sd.2[x])** kernel.power)*weight.2[x]})
  
  
  modal.union.dist.weighted <- lapply(1:length(modal.1.union.dist.weighted), function(x){ modal.1.union.dist.weighted[[x]]+modal.2.union.dist.weighted[[x]] })
  
  
  #select k nearest neighbors - n_neighbors
  select.order <- lapply(modal.union.dist.weighted, function(x){order(x, decreasing = T)}[1:n_neighbors]) #there may be bugs in the original seurat implementation (decreasing = T), seurat implementation is correct
  select.index <- lapply(1:length(select.order), function(x){modal.union.neighbor.index[[x]][select.order[[x]]]})
  select.dist <- lapply(1:length(select.order), function(x){modal.union.dist.weighted[[x]][select.order[[x]]]})
  
  
  weighted.index <- do.call(rbind, select.index)
  weighted.dist <- do.call(rbind, select.dist)
  
  weighted.dist <- sqrt(MinMax((1-weighted.dist)/2, val.min = 0, val.max = 1))
  
  
  #compute KNN
  #########################################
  jj <- c(t(weighted.index))
  ii <- rep(1:nrow(weighted.index), each = ncol(weighted.index))
  
  knn.mat <- Matrix::sparseMatrix(
    i = ii,
    j =jj,
    x = 1,
    dims = c(nrow(weighted.index), nrow(weighted.index))
  )#knn.mat
  
  diag(knn.mat) <- 1
  knn.mat <- knn.mat+Matrix::t(knn.mat)-knn.mat*Matrix::t(knn.mat)
  
  
  
  #compute SNN
  #########################################
  snn.mat <- Seurat:::ComputeSNN(nn_ranked = weighted.index, prune = prune.SNN)
  
  
  
  return(list(weight.1 = weight.1, weight.2 = weight.2, weighted.index = weighted.index, weighted.dist = weighted.dist, knn.mat = knn.mat, snn.mat = snn.mat))
  
}#BiModalIntegration







FindNN <- function(data, query=NULL, number.of.NN, para = BiocNeighbors::AnnoyParam(), prebuild.index = NULL){
  if (!is.null(query)){
    if (is.null(prebuild.index)){
      out <- BiocNeighbors::queryKNN(X = data, query, k = number.of.NN, BNPARAM = para)
    }else{
      out <- BiocNeighbors::queryKNN(query = query, k = number.of.NN, BNPARAM = para, BNINDEX = prebuild.index)  
    }#else
    
  }else{
    if (is.null(prebuild.index)){
      out <- BiocNeighbors::findKNN(X = data, k = number.of.NN, BNPARAM = para)
    }else{
      out <- BiocNeighbors::findKNN(k = number.of.NN, BNPARAM = para, BNINDEX = prebuild.index)
    }#else
  }#else
  
  return(out)
}#FindNN



#assume x is a single vector and y is a matrix with the same number of columns
dist.cal.default <- function(x, y, metric = "euclidean"){
  if (metric=="euclidean"){
    res <- apply(sweep(y,MARGIN = 2, x, "-"), 1, function(x){sqrt(sum(x^2))})
  }#if 
  return(res)
}#dist.cal.default



neighbor_distance_calculation <- function(index.list, mat, metric = "euclidean", nearest.dist = NULL){
  if (nrow(mat)!=length(index.list)){
    stop("Dimension doesn't match.")
  }else{
    neighbor_distance <- lapply(1:length(index.list), function(x){ 
      return( dist.cal.default(mat[x,], mat[index.list[[x]],], metric = metric))
      })#lapply
    
    #print(length(neighbor_distance))
    
    if (!is.null(nearest.dist)){
      if(length(nearest.dist)!=length(index.list)){
        stop("Nearest.dist dimension doesn't match.")
      }else{
      
        neighbor_distance <- lapply(1:length(index.list),function(x){
          
          ndist_adjust <- neighbor_distance[[x]]-nearest.dist[x]
          ndist_adjust <- ReLu(ndist_adjust)
          return(ndist_adjust)
        })#lapply
        
      }#else
    }#if
  }#else
  
  return(neighbor_distance)
}#neighbor_distance_calculation




PredictAssay <- function(mat, neighbor.index){
  mat.pred <- apply(neighbor.index, 1, function(x){ apply(mat[x,,drop=F],2,mean) })
  return(t(mat.pred))
}#PredictAssay



ReLu <- function(x){
  x[x < 0] <- 0
  return(x)
}#ReLu


ImputeDist <- function(x, y, nearest.dist) {
  dist <- sqrt(x = rowSums(x = (x - y)**2)) - nearest.dist
  dist <- ReLu(x = dist)
  return(dist)
}#ImputeDist


MinMax <- function(x, val.min, val.max){
  x[x < val.min] <- val.min
  x[x > val.max] <- val.max
  return(x)
}#MinMax










