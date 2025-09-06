#reference find the neighbor
#https://www.bioconductor.org/packages/release/bioc/vignettes/BiocNeighbors/inst/doc/approx.html





#' Integrate two modalities
#' @import BiocNeighbors
#' @details This function requires `Seurat`. Make sure it is installed
#' @export
BiModalIntegration <- function(
  modal.1,  #to build nearest neighbor 
  modal.2,  #to build nearest neighbor
  mat.1=modal.1,  #the values for prediction
  mat.2=modal.2,  #the values for prediction
  library.size.1 = NULL, 
  library.size.2 = NULL, 
  library.type.1 = "rna", 
  library.type.2 = "atac", 
  baseline.1 = NULL,
  baseline.2 = NULL,
  prior.1 = 1,
  prior.2 = 1,
  n_neighbors = 40, 
  n_neighbors_large = 50, 
  sigma.idx = n_neighbors, 
  snn.far.nn = T, 
  L2norm = "column", 
  sd.scale = 1, 
  cross.contant = 1e-4, 
  prune.SNN = 0, 
  kernel.power = 1){
  
  
  if (sum(duplicated(rownames(modal.1)))>0 | sum(duplicated(rownames(modal.2)))>0){
      stop("The spot/bin names are not unique.")
  }#if
  
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for this function but is not installed. Please install it.")
  }#if

  if ( is.null(rownames(modal.1)) | is.null(rownames(modal.2)) ){
    stop("At least one object has no valid row names. Please examine and assign proper row names to modal.1 and modal.2.")
  }else if ( !all( sort(rownames(modal.1))== sort(rownames(modal.2)) ) ){
    stop("Row names of modal.1 and modal.2 don't align. Please check carefully.")
  }else{
    #make sure modal.1 and modal.2 aligned
    modal.1 <- modal.1[rownames(modal.2),]
    message("Alignment test passed")
  }#else

  message("n_neighbors = ", n_neighbors)
  message("n_neighbors_large = ", n_neighbors_large)
  message("L2 = ", L2norm)
  

  #L2 normalization
  #################################
  if (!is.null(L2norm)){
    if (L2norm=="row"){
    
      message("Sample-wise L2Norm")
      modal.1 <- L2Norm(modal.1, MARGIN = 1)
      modal.2 <- L2Norm(modal.2, MARGIN = 1) 
      mat.1 <- L2Norm(mat.1, MARGIN = 1)
      mat.2 <- L2Norm(mat.2, MARGIN = 1)
      
    }else if (L2norm=="column"){
      
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
  modal.1.neighbor <- FindNN(data = modal.1, number.of.NN = n_neighbors, prebuild.index = prebuild.index.modal.1)
  modal.1.neighbor.large <- FindNN(data = modal.1, number.of.NN = n_neighbors_large, prebuild.index = prebuild.index.modal.1)
  
    
  prebuild.index.modal.2 <- BiocNeighbors::buildIndex(modal.2, BNPARAM = BiocNeighbors::AnnoyParam() )
  modal.2.neighbor <- FindNN(data = modal.2, number.of.NN = n_neighbors, prebuild.index = prebuild.index.modal.2)
  modal.2.neighbor.large <- FindNN(data = modal.2, number.of.NN = n_neighbors_large, prebuild.index = prebuild.index.modal.2)
  
  
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
  #impute.dist.1.NN1 <- ImputeDist(mat.1, mat.1.pred.NN1, NNdist1)
  #impute.dist.2.NN2 <- ImputeDist(mat.2, mat.2.pred.NN2, NNdist2)
  impute.dist.1.NN1 <- ImputeDist(mat.1, mat.1.pred.NN1, 0)
  impute.dist.2.NN2 <- ImputeDist(mat.2, mat.2.pred.NN2, 0)
  #table(impute.dist.1.NN1 > NNdist1)
  #table(impute.dist.2.NN2 > NNdist2)
  
  
  #calculate the distance, across the modality
  #impute.dist.1.NN2 <- ImputeDist(mat.1, mat.1.pred.NN2, NNdist1)
  #impute.dist.2.NN1 <- ImputeDist(mat.2, mat.2.pred.NN1, NNdist2)
  impute.dist.1.NN2 <- ImputeDist(mat.1, mat.1.pred.NN2, 0)
  impute.dist.2.NN1 <- ImputeDist(mat.2, mat.2.pred.NN1, 0)
  #table(impute.dist.1.NN2 > NNdist1)
  #table(impute.dist.2.NN1 > NNdist2)
  
  
  #calculate the kernel width
  #################################
  if (snn.far.nn){
    snn.matrix.1 <- Seurat:::ComputeSNN(nn_ranked = modal.1.neighbor$index, prune = prune.SNN)
    snn.matrix.2 <- Seurat:::ComputeSNN(nn_ranked = modal.2.neighbor$index, prune = prune.SNN)
    
    #snn.width.1 <- Seurat:::ComputeSNNwidth(snn.graph = snn.matrix.1, k.nn = n_neighbors, l2.norm = F, embeddings = mat.1, nearest.dist = NNdist1)
    
    snn.width.1 <- Seurat:::ComputeSNNwidth(
      snn.graph = snn.matrix.1, 
      k.nn = n_neighbors, 
      l2.norm = F, 
      embeddings = mat.1, 
      nearest.dist = rep(0, nrow(modal.1))
      )
    
    #snn.width.2 <- Seurat:::ComputeSNNwidth(snn.graph = snn.matrix.2, k.nn = n_neighbors, l2.norm = F, embeddings = mat.2, nearest.dist = NNdist2)
    snn.width.2 <- Seurat:::ComputeSNNwidth(
      snn.graph = snn.matrix.2, 
      k.nn = n_neighbors, 
      l2.norm = F, 
      embeddings = mat.2, 
      nearest.dist = rep(0, nrow(modal.2))
      )
    
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
  
  #score.1 <- 2*score.1
  #score.1 <- 2*score.2
  
  
  score.sum <- exp(score.1)+exp(score.2)
  weight.1 <- exp(score.1)/score.sum
  weight.2 <- exp(score.2)/score.sum
  
  

  #reweigh the weight with the library size
  #################################
  #library.size.1 will be a vector
  #library.size.2 will be a vector
  message("")
  if ( !is.null(library.size.1) & !is.null(library.size.2) ){
      
       print("Adjust the weights by library size.")
       
       #with internel parameters
       #calculated based on /mnt/ceph/users/wmao1/KPMP/code/PLIER_Ahmed/speed_sylvester/exhaustive/CoProfiling/MouseBrain/generate_input.R
       #rna per spot
       baseline.rna <- 2300 #mean, can think of median
       #atac per spot
       baseline.atac <- 7400

      if (is.null(baseline.1)){

          if (library.type.1 == "rna"){
                baseline.1 <- baseline.rna
          }else if (library.type.1 == "atac"){
                baseline.1 <- baseline.atac
          }#else if

          message("Reference UMI of modality 1 is set to be: ", baseline.1)
      }#if is.null
      
      if (is.null(baseline.2)){
              if (library.type.2 == "rna"){
                      baseline.2 <- baseline.rna
              }else if (library.type.2 == "atac"){
                      baseline.2 <- baseline.atac
              }#else if 

              message("Reference UMI of modality 2 is set to be: ", baseline.2)
      }#if is.null

      
      #with thresholding
      library.ratio.1 <- pmin(library.size.1/baseline.1, 1)
      library.ratio.2 <- pmin(library.size.2/baseline.2, 1)
      
      
      #quantile thresholding to get more dramatic results
      #library.ratio.1 <- ecdf(library.ratio.1)(library.ratio.1)
      #library.ratio.2 <- ecdf(library.ratio.2)(library.ratio.2)
      
      
      #deno <- exp(weight.1*library.ratio.1)+exp(weight.2*library.ratio.2)
      #weight.1 <- exp(weight.1*library.ratio.1)/deno
      #weight.2 <- exp(weight.2*library.ratio.2)/deno

      weight.1 <- weight.1*library.ratio.1
      weight.2 <- weight.2*library.ratio.2
      
      deno <- weight.1+weight.2
      weight.1 <- weight.1/deno
      weight.2 <- weight.2/deno

  }else{
          print("No adjustment with reference UMI.")
  }#else


#apply the prior re-weighting - numerical stability
#weight.1 <- weight.1*prior.1
#weight.2 <- weight.2*prior.2
#deno <- weight.1+weight.2
#weight.1 <- weight.1/deno
#weight.2 <- weight.2/deno




  message("Calculating the weighted KNN and SNN")
  ########################################################################################
  
  if (T){
    
    #calculate the weighted sum of the link weights
    n_neighbors_k <- n_neighbors
 
    modal.union.neighbor.index <- list()
    for (ii in 1:nrow(modal.1.neighbor.large$index)){
      
      index.intersection <- intersect(modal.1.neighbor.large$index[ii,], modal.2.neighbor.large$index[ii,])
      index.1.only <- setdiff(modal.1.neighbor.large$index[ii,], index.intersection)
      index.2.only <- setdiff(modal.2.neighbor.large$index[ii,], index.intersection)
      
      index.union <- c(index.intersection, index.1.only, index.2.only)
      index.union.lw <- c(weight.1[index.intersection]+weight.2[index.intersection],
                          weight.1[index.1.only],
                          weight.2[index.2.only]
                          )
      
      #sort and take the maxmimal several
      modal.union.neighbor.index[[ii]] <- index.union[order(index.union.lw, decreasing = T)[1:n_neighbors_k]]
        
    }#for ii
 
    weighted.index <- do.call(rbind, modal.union.neighbor.index)
    weighted.dist <- modal.union.neighbor.index
    
    
    
  }else{
    
    #figure out the union set of neighbors - n_neighbors_large
    modal.union.neighbor.index <- list()
    for (ii in 1:nrow(modal.1.neighbor.large$index)){
      modal.union.neighbor.index[[ii]] <- union(modal.1.neighbor.large$index[ii,], modal.2.neighbor.large$index[ii,])
    }#for ii
    
    
    
    #calculate the new distance matrix in each modality given the union set
    #modal.1.union.dist <- neighbor_distance_calculation(modal.union.neighbor.index, modal.1, nearest.dist = NNdist1)
    #modal.2.union.dist <- neighbor_distance_calculation(modal.union.neighbor.index, modal.2, nearest.dist = NNdist2)
    modal.1.union.dist <- neighbor_distance_calculation(modal.union.neighbor.index, modal.1, nearest.dist = rep(0, nrow(modal.1)))
    modal.2.union.dist <- neighbor_distance_calculation(modal.union.neighbor.index, modal.2, nearest.dist = rep(0, nrow(modal.2)))
    
    
    
    #weighted by the modal weights, and then take the total sum
    modal.1.union.dist.weighted <- lapply(1:length(modal.1.union.dist), function(x){ exp(-1*(modal.1.union.dist[[x]]/sd.1[x])** kernel.power)*weight.1[x]})
    modal.2.union.dist.weighted <- lapply(1:length(modal.2.union.dist), function(x){ exp(-1*(modal.2.union.dist[[x]]/sd.2[x])** kernel.power)*weight.2[x]})
    modal.union.dist.weighted <- lapply(1:length(modal.1.union.dist.weighted), function(x){ modal.1.union.dist.weighted[[x]]+modal.2.union.dist.weighted[[x]] })
    
    
    
    #select k nearest neighbors - n_neighbors
    n_neighbors_k <- min(unlist(lapply(modal.union.neighbor.index, length)))
    print(n_neighbors_k)
    if (n_neighbors_k < n_neighbors){
      n_neighbors_k <- n_neighbors
    }#if
    
    select.order <- lapply(modal.union.dist.weighted, function(x){order(x, decreasing = T)}[1:n_neighbors_k]) #there may be bugs in the original seurat implementation (decreasing = T), - seurat implementation is correct
    select.index <- lapply(1:length(select.order), function(x){modal.union.neighbor.index[[x]][select.order[[x]]]})
    select.dist <- lapply(1:length(select.order), function(x){modal.union.dist.weighted[[x]][select.order[[x]]]})
    
    
    weighted.index <- do.call(rbind, select.index)
    
    weighted.dist <- do.call(rbind, select.dist)
    weighted.dist <- sqrt(MinMax((1-weighted.dist)/2, val.min = 0, val.max = 1))
    
  }#else
  
  
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
  
  
  
  return(list(weight.1 = weight.1, 
              weight.2 = weight.2, 
              weighted.index = weighted.index, 
              weighted.dist = weighted.dist, 
              knn.mat = knn.mat, 
              snn.mat = snn.mat,
              library.size.1 = library.size.1, 
              library.size.2 = library.size.2, 
              library.type.1 = library.type.1, 
              library.type.2 = library.type.2, 
              baseline.1 = baseline.1,
              baseline.2 = baseline.2,
              prior.1 = prior.1,
              prior.2 = prior.2,
              n_neighbors = n_neighbors, 
              n_neighbors_large = n_neighbors_large))
  
}#BiModalIntegration
















BiModalIntegration.legacy <- function(modal.1, 
                              modal.2, 
                              mat.1=modal.1, 
                              mat.2=modal.2, 
                              n_neighbors = 20, 
                              n_neighbors_large = 200, 
                              sigma.idx = n_neighbors, 
                              snn.far.nn = T, 
                              L2norm = "col", 
                              sd.scale = 1, 
                              cross.contant = 1e-4, 
                              prune.SNN = 0, 
                              kernel.power = 1){
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for this function but is not installed. Please install it.")
  }#if

  if ( is.null(rownames(modal.1)) | is.null(rownames(modal.2)) ){
    stop("At least one object has no valid row names. Please examine and assign proper row names to modal.1 and modal.2.")
  }else if ( !all( sort(rownames(modal.1))== sort(rownames(modal.2)) ) ){
    stop("Row names of modal.1 and modal.2 don't align. Please check carefully.")
  }else{
    #make sure modal.1 and modal.2 aligned
    modal.1 <- modal.1[rownames(modal.2),]
    message("Alignment test passed")
  }#else


  #L2 normalization
  #################################
  if (!is.null(L2norm)){
    if (L2norm == "row"){
    
      message("Sample-wise L2Norm")
      modal.1 <- L2Norm(modal.1, MARGIN = 1)
      modal.2 <- L2Norm(modal.2, MARGIN = 1) 
      mat.1 <- L2Norm(mat.1, MARGIN = 1)
      mat.2 <- L2Norm(mat.2, MARGIN = 1)
      
    }else if (L2norm == "column"){
      
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
  
  

  #reweigh the weight with the library size
  #################################

  #with thresholding

  
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
  
}#BiModalIntegration.legacy







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










