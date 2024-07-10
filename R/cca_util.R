
impact.default <- function(reference.Z, query.gene.exp, k, L2, right.shur, thr = 0.8){
  
  left <- 2*crossprod(reference.Z)+2*L2*diag(k)
  total <- 2*t(reference.Z) %*% query.gene.exp
  
  B <- sylvester_pre(left, right.shur$U, right.shur$S, total)
  
  #L4_adaptive
  ###########################################
  B.cor.res <- cor(t(B))
  B.cor.val <- B.cor.res[upper.tri(B.cor.res)]
  #pheatmap::pheatmap(B.cor.res, display_numbers = T, fontsize = 15)
  
  
  #print(sum(B.cor.val>=thr))
  if (sum(B.cor.val>=thr)==1){
    flag <- "Done"
  }else if (sum(B.cor.val>=thr)>1){
    flag <- "Shrink"
  }else{
    flag <- "Explode"
  }#else
  
  
  return(list(flag = flag, B = B))  
}#impact.default





#update L4
impact_adaptive <- function(reference.Z, query.gene.exp, query.L, query.L1 = NULL, query.L2 = NULL, query.L4, query.shur0 = NULL, query.ICAp.res = NULL, to_drop = T, scale=1, L4_adaptive =2, max.iter = 200, cor.thr = 0.8, save.complete=T, verbose = T){
  pos.adj <- 3
  
  k <- ncol(reference.Z)
  
  #no need
  if (to_drop){
    k <- k+1
  }#if
  
  #align the reference.Z and query.gene.exp
  gene.intersect <- intersect(rownames(reference.Z), rownames(query.gene.exp))
  reference.Z <- reference.Z[gene.intersect,]
  query.gene.exp <- query.gene.exp[gene.intersect,]
  
  
  if (verbose){
    message("************")  
  }#if verbose
  
  
  #calculate the L1, L2 parameters
  ###########################################################
  if (is.null(query.ICAp.res)){
    if (is.null(query.L1) | is.null(query.L2)){
      
      if (verbose){
        message("Computing SVD")  
      }#if verbose
      
      set.seed(123)
      svdres <- rsvd(query.gene.exp, k = k)
      svdres <- rotateSVD(svdres)
      
      L1 <- svdres$d[k]*scale
      if (!is.null(pos.adj)){
        L1 <- L1/pos.adj
      }#if
      L2 <- svdres$d[k]*scale
      
    }else{
    
      L1 <- query.L1
      L2 <- query.L2
      
    }#else
    
  }else{
    L1 <- query.ICAp.res$L1
    L2 <- query.ICAp.res$L2
  }#else
  
  if (verbose){
    print(paste0("L1 is set to ", L1))
    print(paste0("L2 is set to ", L2))  
  }#if verbose
  
  
  L4 <- query.L4
  #B
  #######################################################################
  
  if (is.null(query.shur0)){
    query.shur0 <- rcpp_shur(query.L)  
  }#if
  
  right.shur <- query.shur0
  right.shur$S <- 2*L4*right.shur$S
  
  if (verbose){
    message("Shur done")  
  }#if verbose
  
  md.run <- impact.default(reference.Z, query.gene.exp, k, L2, right.shur, thr = cor.thr)
  
  if (md.run$flag == "Explode"){
    if (verbose){
      message("Explode")  
    }#if verbose
    
    
    L4_left <- L4
    L4_right <- L4*L4_adaptive
    L4_pointer <- L4_right
  }else if (md.run$flag =="Shrink"){
    
    if (verbose){
      message("Shrink")  
    }#if verbose
    
    L4_left <- L4/L4_adaptive
    L4_right <- L4
    L4_pointer <- L4_left
  }else{
    L4_pointer <- L4
    L4_left <- NULL
    L4_right <- NULL
  }#else
  
  if (verbose){
    cat(paste0(L4_pointer, " ", L4_left, " ", L4_right, "\n"))  
  }#if verbose
  
  
  while (md.run$flag != "Done"){
    
    #right <- 2*L4_pointer*L
    #right.shur <- rcpp_shur(right)
    right.shur <- query.shur0
    right.shur$S <- 2*L4_pointer*right.shur$S
    
    #if L4_left and L4_right is too close, then Done
    if (abs(L4_left-L4_right) < 0.1){
      md.run <- impact.default(reference.Z, query.gene.exp, k, L2, right.shur, thr = cor.thr)
      md.run$flag <- "Done"
    }else{
      md.run <- impact.default(reference.Z, query.gene.exp, k, L2, right.shur, thr = cor.thr)
    }#else
    
    
    if (md.run$flag == "Explode"){
      if (verbose){
        message("Explode")  
      }#if verbose
      
      if (L4_pointer == L4_left){
        L4_left <- L4_left 
        L4_pointer <- (L4_left+L4_right)/2
        L4_right <- L4_right
        
      }else if (L4_pointer < L4_right){
        L4_left <- L4_pointer
        L4_pointer <- (L4_pointer+L4_right)/2
        L4_right <- L4_right
        
      }else if (L4_pointer == L4_right){
        L4_pointer <- L4_right*L4_adaptive
        L4_left <- L4_right
        L4_right <- L4_pointer
      }#else if 
      
    }else if (md.run$flag == "Shrink"){
      if (verbose){
        message("Shrink")  
      }#if verbose
      
      if (L4_pointer == L4_left){
        L4_pointer <- L4_pointer/L4_adaptive
        L4_right <- L4_left
        L4_left <- L4_pointer
        
      }else if (L4_pointer < L4_right){
        L4_right <- L4_pointer
        L4_pointer <- (L4_left+L4_pointer)/2
        L4_left <- L4_left
      }else if (L4_pointer == L4_right){
        L4_left <- L4_left
        L4_pointer <- (L4_left+L4_right)/2
        L4_right <- L4_right
      }#else if
      
    }#else if 
    
    if (verbose){
      cat(paste0(L4_pointer, " ", L4_left, " ", L4_right, "\n"))  
    }#if verbose
    
    
  }#while
  
  
  #wrap around the output
  ################################################################
  B <- md.run$B
  
  if (to_drop){
    
    if (verbose){
      message("drop")  
    }#if verbose
    
    #eliminate the one with smaller variance
    cor.res <- cor(t(B))
    LV.var <- VarianceExplained(query.gene.exp, reference.Z,B, option="simple", normalize = F)
    drop.index <- which(cor.res == max(cor.res[upper.tri(cor.res)]) & upper.tri(cor.res), arr.ind = T)
    if (nrow(drop.index)!=1){
      warning(nrow(drop.index))
    }#if
    LV.to_drop <- drop.index[which.min(LV.var[drop.index])]
    
    #exclue LV.to_drop
    B <- B[-LV.to_drop,]
    k <- k-1
  }else{
    LV.to_drop <- NULL
  }#else
  
  rownames(B) <- paste0("LV",1:k)
  colnames(B) <- colnames(query.gene.exp)
  
  if (save.complete){
    return(list(B=B, L1 = L1, L2 = L2, L4 = L4_pointer, k = k, shur0 = query.shur0, right.shur = right.shur, reference.Z = reference.Z, query.gene.exp = query.gene.exp, query.L = query.L, LV.to_drop = LV.to_drop))
  }else{
    return(list(B=B, L1 = L1, L2 = L2, L4 = L4_pointer, k = k, shur0 = query.shur0, right.shur = right.shur))
  }#else
  
}#impact_adaptive




#use the fixed L4
impact <- function(reference.Z, query.gene.exp, query.L, query.L4, query.shur0 = NULL, query.ICAp.res = NULL, scale=1, max.iter = 200, cor.thr = 0.8, save.complete=T){
  pos.adj <- 3
  
  k <- ncol(reference.Z)
  
  #align the reference.Z and query.gene.exp
  gene.intersect <- intersect(rownames(reference.Z), rownames(query.gene.exp))
  reference.Z <- reference.Z[gene.intersect,]
  query.gene.exp <- query.gene.exp[gene.intersect,]
  
  
  message("************")
  
  #calculate the L1, L2 parameters
  ###########################################################
  if (is.null(query.ICAp.res)){
    message("Computing SVD")
    set.seed(123)
    svdres <- rsvd(query.gene.exp, k = k)
    svdres <- rotateSVD(svdres)
    
    L1 <- svdres$d[k]*scale
    if (!is.null(pos.adj)){
      L1 <- L1/pos.adj
    }#if
    L2 <- svdres$d[k]*scale
    
  }else{
    L1 <- query.ICAp.res$L1
    L2 <- query.ICAp.res$L2
  }#else
  
  print(paste0("L1 is set to ", L1))
  print(paste0("L2 is set to ", L2))
  
  L4 <- query.L4
  #B
  #######################################################################
  
  if (is.null(query.shur0)){
    query.shur0 <- rcpp_shur(query.L)  
  }#if
  
  right.shur <- query.shur0
  right.shur$S <- 2*L4*right.shur$S
  message("Shur done")
  
  md.run <- impact.default(reference.Z, query.gene.exp, k, L2, right.shur, thr = cor.thr)
  
  
  #wrap around the output
  ################################################################
  B <- md.run$B
  
  rownames(B) <- paste0("LV",1:k)
  colnames(B) <- colnames(query.gene.exp)
  
  if (save.complete){
    return(list(B=B, L1 = L1, L2 = L2, L4 = L4, k = k, query.shur0 = query.shur0, right.shur = right.shur, reference.Z = reference.Z, query.gene.exp = query.gene.exp, query.L = query.L))
  }else{
    return(list(B=B, L1 = L1, L2 = L2, L4 = L4, k = k, query.shur0 = query.shur0, right.shur = right.shur))
  }#else
  
}#impact







Z.align <- function(Z1, Z2, plot = T){
  align.res <- cor.align(Z1, Z2, B.list=NULL)
  
  gene.intersect <- intersect(rownames(Z1), rownames(Z2))
  Z1 <- Z1[gene.intersect,]
  Z2 <- Z2[gene.intersect,]
    
  cor.res <- cor(Z1, Z2)

  if (plot){
  pheatmap::pheatmap(cor.res[,align.res$adjust.index], display_numbers = T, fontsize = 20, cluster_rows = F, cluster_cols = F)
  }#if plot
  
  return(list(cor.res = cor.res[,align.res$adjust.index] , adjust.index = align.res$adjust.index))
}#Z.align





cor.align <- function(mat1, mat2, B.list=NULL){

  gene.intersect <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[gene.intersect,]
  mat2 <- mat2[gene.intersect,]
  
  
  cor.res <- cor(mat1, mat2) 
  cor.res <- abs(cor.res)
  
  #fix 1 and find order in 2
  mat1.index <- c()
  mat2.index <- c()
  
  mat1.index.available <- 1:ncol(mat1)
  mat2.index.available <- 1:ncol(mat2)
  
  while (length(mat1.index.available) >0){
    temp <- which(cor.res == max(cor.res), arr.ind = T)
    mat1.index <- c(mat1.index, temp[1,1])
    mat2.index <- c(mat2.index, temp[1,2])
    
    #adjust cor.res and other two variables
    mat1.index.available <- setdiff(mat1.index.available, mat1.index)
    mat2.index.available <- setdiff(mat2.index.available, mat2.index)
    cor.res[temp[1,1],] <- -1
    cor.res[,temp[1,2]] <- -1
  }#while

  
  B.adjust <-  B.list
  if (!is.null(B.list)){
    B.adjust[[2]] <- B.adjust[[2]][mat2.index[order(mat1.index)],]
  }#if
  
  return(list(B.adjust = B.adjust, adjust.index = mat2.index[order(mat1.index)]))
  
}#cor.align



#unfortunately, the implemenration here (#reference: /mnt/ceph/users/wmao1/KPMP/code/mirror/mirror.R) is wrong
mirror <- function(Y, B, Lap, Pairing, L2, L3, k, max.iter=200, tol=5e-6, trace = F){
  #Y, B, Lap are lists
  #L2, L3 are vectors
  #Pairing is a list of vectors
  
  
  first <- Pairing[[1]][1]
  second <- Pairing[[1]][2]
  
  
  ptm <- proc.time()
  #left_1 <- t(Y[[first]])%*%Y[[first]]+L2[first]*diag(ncol(Y[[first]]))+L3[first]*Lap[[first]]
  left_1 <- eigenMapMatMult(t(Y[[first]]), Y[[first]])+L2[first]*diag(ncol(Y[[first]]))+L3[first]*Lap[[first]]
  left_1 <- solve(left_1)
  #right_1 <- t(Y[[second]])%*% Y[[first]]
  right_1 <- eigenMapMatMult(t(Y[[second]]), Y[[first]])
  
  
  #left_2 <- t(Y[[second]])%*%Y[[second]]+L2[second]*diag(ncol(Y[[second]]))+L3[second]*Lap[[second]]
  left_2 <- eigenMapMatMult(t(Y[[second]]), Y[[second]])+L2[second]*diag(ncol(Y[[second]]))+L3[second]*Lap[[second]]
  left_2 <- solve(left_2)
  #right_2 <- t(Y[[first]]) %*% Y[[second]]
  right_2 <- eigenMapMatMult(t(Y[[first]]), Y[[second]])
  print(proc.time()-ptm)
  
  
  BdiffTrace <- c()
  BdiffCount <- 0
  
  for (ii in 1:max.iter){
    
    oldB <- B
    
    #update each B in order
    ####################################################
    #TBD: multiple
    for (jj in 1:length(Pairing)){
      first <- Pairing[[jj]][1]
      second <- Pairing[[jj]][2]
      
      #Update B_first
      #left <- t(Y[[first]])%*%Y[[first]]+L2[first]*diag(ncol(Y[[first]]))+L3[first]*Lap[[first]]
      #right <- B[[second]]%*% t(Y[[second]])%*% Y[[first]]
      right <- B[[second]]%*%right_1
      B[[first]] <- right%*% left_1
      #print("first")
      
      #Update B_second
      #left <- t(Y[[second]])%*%Y[[second]]+L2[second]*diag(ncol(Y[[second]]))+L3[second]*Lap[[second]]
      #right <- B[[first]] %*% t(Y[[first]]) %*% Y[[second]]
      right <- B[[first]] %*%right_2
      B[[second]] <- right%*%left_2
      #print("second")  
        
    }#for jj
    
    #update error
    ####################################################
    #calculate the difference
    Bdiff <- c()
    for (jj in 1:length(B)){
      Bdiff <- c(Bdiff, sum((B[[jj]]-oldB[[jj]])^2)/sum(B[[jj]]^2))
    }#for jj
    
    Bdiff.mean <- mean(Bdiff)
    BdiffTrace <- c(BdiffTrace, Bdiff.mean)
    
    #calculate the error
    ########################
    if (trace & (ii%%1==0) ){
      print(error(Y, B, Pairing, Lap, L2, L3))
      print(Bdiff.mean)
      #print(BdiffTrace)
    }#if 
    
    #check the convergence
    ########################
    if(ii>52){
      if (Bdiff.mean > BdiffTrace[ii-50]){
        BdiffCount <- BdiffCount+1  
      }else if(BdiffCount>1){
        BdiffCount <- BdiffCount-1
      }#else if 
    }else if(BdiffCount>1){
      BdiffCount <- BdiffCount-1
    }#else if 
    
    #if(Bdiff.mean < tol & ii>40){
    if(Bdiff.mean < tol & ii>20){
      message(paste0("converged at  iteration ", ii))
      break
    }#if
    
    #if( BdiffCount>5 & ii>40){
    if( BdiffCount>5 & ii>20){
      message(paste0("stopped at  iteration ", ii, " Bdiff.mean is not decreasing"))
      break
    }#if
    
  }#for ii
  
  return(list(B=B, L2=L2, L3 = L3, k = k, number.of.iter = ii, BdiffTrace = BdiffTrace))
}#mirror




tr <- function(mat){
  return(sum(diag(mat)))
}#tr



error <- function(Y, B, Pairing, Lap, L2, L3){
  nn <- length(Y)
  res <- 0
  
  for (ii in 1:length(Pairing)){
    first <- Pairing[[ii]][1]
    second <- Pairing[[ii]][2]
    res <- res+normF( Y[[first]]%*%t(B[[first]])- Y[[second]]%*%t(B[[second]]) )
  }#for ii
  
  for (ii in 1:nn){
    res <- res+L2[[ii]]*normF(B[[ii]])+L3[[ii]]*tr(B[[ii]] %*% Lap[[ii]] %in% t(B[[ii]]) )
  }#for i
  
  return(res)
}#error






parse <- function(input_path, slice_name){
  dat <- readRDS(paste0(input_path, slice_name, ".RDS"))
  gene.exp <- dat@assays$Spatial@counts
  
  
  use_gene <- which(apply(gene.exp==0,1,sum)/ncol(gene.exp) < 0.95)
  gene.exp <- as.matrix(gene.exp[use_gene,])
  
  
  #remove MT- genes
  use_gene <- which(!grepl("MT-", rownames(gene.exp)))
  gene.exp <- gene.exp[use_gene,]
  
  
  gene.exp <- normalize_barcode_sums_to_median(gene.exp)
  gene.exp <- log2(gene.exp+1)
  
  #prepare the laplacian
  if (all( c("array_row", "array_col") %in% colnames(dat@meta.data))){
    coor <- dat@meta.data[,c("array_row", "array_col")]  
  }else{
    coor <- dat@images$slice1@coordinates[,c("row", "col")]  
  }#else
  
  coor.dist <- dist(coor)
  coor.dist <- as.matrix(coor.dist)
  
  
  x <- c()
  y <- c()
  x.end <- c()
  y.end <- c()
  A <- matrix(0, nrow = nrow(coor), ncol = nrow(coor))
  rownames(A) <- rownames(coor)
  colnames(A) <- rownames(coor)
  
  
  for (i in 1: (nrow(coor.dist)-1)){
    temp <- which(coor.dist[i, (i+1):ncol(coor.dist)]<=2)
    if (length(temp)>0){
      A[rownames(coor.dist)[i], names(temp)] <- 1
      A[names(temp), rownames(coor.dist)[i]] <- 1
      x <- c(x, rep(coor[i,1], length(temp)))
      y <- c(y, rep(coor[i,2], length(temp)))
      x.end <- c(x.end, coor[names(temp),1])
      y.end <- c(y.end, coor[names(temp),2])
    }#if
  }#for i
  
  D <-  diag(apply(A,1,sum))
  L <- D-A
  
  return(list(gene.exp=gene.exp, L=L, dat = dat))
}#parse






prepare <- function(input.list){
  gene.intersect <- rownames(input.list[[1]]$gene.exp)
  
  for (ii in 2:length(input.list)){
    gene.intersect <- intersect(gene.intersect, rownames(input.list[[ii]]$gene.exp))    
  }#for ii
  
  Y <- list()
  B <- list()
  Lap <- list()
  
  for (ii in 1:length(input.list)){
    Y[[ii]] <- input.list[[ii]]$gene.exp[gene.intersect,]
    Lap[[ii]] <- input.list[[ii]]$L
    B[[ii]] <- input.list[[ii]]$ICAp.res$B
  }#for ii
  
  
  count <- 0
  Pairing <- list()
  for (ii in 1:(length(input.list)-1)){
    for (jj in (ii+1):length(input.list)){
      count <- count+1
      Pairing[[count]] <- c(ii, jj)
    }#for jj
  }#for ii
  
  return(list(Y = Y, Lap = Lap, B = B, Paring = Pairing))
}#prepare













constrained_mirror <- function(Y.list, Lap, ICAp.res.list, right.shur.list, L3.list, k, max.iter=200, tol=5e-6, trace = F){
  #Y, Lap and ICAp.res.list are lists
  
  #align the genes
  common.genes <- intersect(rownames(Y.list[[1]]), rownames(Y.list[[2]]))
  Y.list[[1]] <- Y.list[[1]][common.genes,]
  Y.list[[2]] <- Y.list[[2]][common.genes,]
  
  ICAp.res.list[[1]]$Z <- ICAp.res.list[[1]]$Z[common.genes, ]
  ICAp.res.list[[2]]$Z <- ICAp.res.list[[2]]$Z[common.genes, ]
  
  k <- nrow(ICAp.res.list[[1]]$B)
  
  
  BdiffTrace <- c()
  BdiffCount <- 0
  
  
  for (ii in 1:max.iter){
  
    oldB <- list(ICAp.res.list[[1]]$B, ICAp.res.list[[2]]$B)
    
    ######################################################################################################################
    for (jj in 1:length(Y.list)){
      Y <- Y.list[[jj]]
      Z <- ICAp.res.list[[jj]]$Z
      B <- ICAp.res.list[[jj]]$B
      Z0 <- ICAp.res.list[[3-jj]]$Z
      L <- Lap[[jj]]
      L1 <- ICAp.res.list[[jj]]$L1
      L2 <- ICAp.res.list[[jj]]$L2
      L3 <- L3.list[[jj]]
      L4 <- ICAp.res.list[[jj]]$L4
      right.shur <- right.shur.list[[jj]]
      
      #Z update
      ##################
      Zraw=Z=(eigenMapMatMult(Y, t(B))+L3*Z0)%*%solve(tcrossprod(B)+(L1+L3)*diag(k))
      Z[Z<0]=0
      
      
      #B update
      ##################
      left <- 2*crossprod(Z)+2*L2*diag(k)
      total <- 2*t(Z) %*% Y
      B <- sylvester_pre(left, right.shur$U, right.shur$S, total)
      
      
      #assign the value back
      ##################
      ICAp.res.list[[jj]]$Z <- Z
      ICAp.res.list[[jj]]$Zraw <- Zraw
      ICAp.res.list[[jj]]$B <- B
      
    }#for jj
    ######################################################################################################################
    
    
    #update error
    ####################################################
    #calculate the difference
    Bdiff <- c()
    newB <- list(ICAp.res.list[[1]]$B, ICAp.res.list[[2]]$B)
    for (jj in 1:length(newB)){
      Bdiff <- c(Bdiff, sum((newB[[jj]]-oldB[[jj]])^2)/sum(oldB[[jj]]^2))
    }#for jj
    
    Bdiff.mean <- mean(Bdiff)
    BdiffTrace <- c(BdiffTrace, Bdiff.mean)
    
    #calculate the error
    ########################
    if (trace & (ii%%1==0) ){
      print(Bdiff.mean)
      #print(BdiffTrace)
    }#if 
    
    #check the convergence
    ########################
    if(ii>52){
      if (Bdiff.mean > BdiffTrace[ii-50]){
        BdiffCount <- BdiffCount+1  
      }else if(BdiffCount>1){
        BdiffCount <- BdiffCount-1
      }#else if 
    }else if(BdiffCount>1){
      BdiffCount <- BdiffCount-1
    }#else if 
    
    #if(Bdiff.mean < tol & ii>40){
    if(Bdiff.mean < tol & ii>20){
      message(paste0("converged at  iteration ", ii))
      break
    }#if
    
    #if( BdiffCount>5 & ii>40){
    if( BdiffCount>5 & ii>20){
      message(paste0("stopped at  iteration ", ii, " Bdiff.mean is not decreasing"))
      break
    }#if
    
  }#for ii
  
  #assign names
  for (jj in 1:length(Y.list)){
    rownames(newB[[jj]]) = rownames(ICAp.res.list[[jj]]$B) = paste0("LV",1:k)
    colnames(newB[[jj]]) = colnames(ICAp.res.list[[jj]]$B) = colnames(Y.list[[jj]])
    rownames(ICAp.res.list[[jj]]$Z) = rownames(Y.list[[jj]])
    colnames(ICAp.res.list[[jj]]$Z) = paste0("LV",1:k)
  }#for jj
  
  return(list(ICAp.res.list=ICAp.res.list, B=newB, k = k, number.of.iter = ii, BdiffTrace = BdiffTrace))
}#constrained_mirror






gene.align <- function(dat.list){
  common.genes <- rownames(dat.list[[1]])
  for (i in 2:length(dat.list)){
    common.genes <- intersect(common.genes, rownames(dat.list[[i]]))
  }#for i
  return(common.genes)
}#gene.align



multi_mirror <- function(Y.list, Lap, ICAp.res.list, right.shur.list, L3.list, k, max.iter=200, tol=5e-6, trace = F){
  #Y, Lap and ICAp.res.list are lists
  k <- nrow(ICAp.res.list[[1]]$B)
  number.of.slice <- length(Y.list)
  
  
  #align the genes
  common.genes <- gene.align(Y.list)
  for (ii in 1:number.of.slice){
    Y.list[[ii]] <- Y.list[[ii]][common.genes,]
    ICAp.res.list[[ii]]$Z <- ICAp.res.list[[ii]]$Z[common.genes, ]
  }#for ii
  
  
  
  BdiffTrace <- c()
  BdiffCount <- 0
  
  
  for (ii in 1:max.iter){
    
    oldB <- lapply(ICAp.res.list, function(x){x$B})
    Z.total <- do.call(sum, lapply(ICAp.res.list, function(x){x$Z}))
    
    ######################################################################################################################
    for (jj in 1:length(Y.list)){
      Y <- Y.list[[jj]]
      Z <- ICAp.res.list[[jj]]$Z
      B <- ICAp.res.list[[jj]]$B
      Z0 <- Z.total-ICAp.res.list[[jj]]$Z
      L <- Lap[[jj]]
      L1 <- ICAp.res.list[[jj]]$L1
      L2 <- ICAp.res.list[[jj]]$L2
      L3 <- L3.list[[jj]]
      L4 <- ICAp.res.list[[jj]]$L4
      right.shur <- right.shur.list[[jj]] #already take L4 into consideration
      
      #Z update
      ##################
      Zraw=Z=(eigenMapMatMult(Y, t(B))+L3*Z0)%*%solve(tcrossprod(B)+(L1+L3*(number.of.slice-1))*diag(k))
      Z[Z<0]=0
      
      
      #B update
      ##################
      left <- 2*crossprod(Z)+2*L2*diag(k)
      total <- 2*t(Z) %*% Y
      B <- sylvester_pre(left, right.shur$U, right.shur$S, total)
      
      
      #assign the value back
      ##################
      ICAp.res.list[[jj]]$Z <- Z
      ICAp.res.list[[jj]]$Zraw <- Zraw
      ICAp.res.list[[jj]]$B <- B
      
      #update Z.total
      ##################
      Z.total <- Z0+Z
      
    }#for jj
    ######################################################################################################################
    
    
    #update error
    ####################################################
    #calculate the difference
    Bdiff <- c()
    newB <- lapply(ICAp.res.list, function(x){x$B})
    
    for (jj in 1:length(newB)){
      Bdiff <- c(Bdiff, sum((newB[[jj]]-oldB[[jj]])^2)/sum(oldB[[jj]]^2))
    }#for jj
    
    Bdiff.mean <- mean(Bdiff)
    BdiffTrace <- c(BdiffTrace, Bdiff.mean)
    
    #calculate the error
    ########################
    if (trace & (ii%%1==0) ){
      print(Bdiff.mean)
      #print(BdiffTrace)
    }#if 
    
    #check the convergence
    ########################
    if(ii>52){
      if (Bdiff.mean > BdiffTrace[ii-50]){
        BdiffCount <- BdiffCount+1  
      }else if(BdiffCount>1){
        BdiffCount <- BdiffCount-1
      }#else if 
    }else if(BdiffCount>1){
      BdiffCount <- BdiffCount-1
    }#else if 
    
    #if(Bdiff.mean < tol & ii>40){
    if(Bdiff.mean < tol & ii>20){
      message(paste0("converged at  iteration ", ii))
      break
    }#if
    
    #if( BdiffCount>5 & ii>40){
    if( BdiffCount>5 & ii>20){
      message(paste0("stopped at  iteration ", ii, " Bdiff.mean is not decreasing"))
      break
    }#if
    
  }#for ii
  
  #assign names
  for (jj in 1:length(Y.list)){
    rownames(newB[[jj]]) = rownames(ICAp.res.list[[jj]]$B) = paste0("LV",1:k)
    colnames(newB[[jj]]) = colnames(ICAp.res.list[[jj]]$B) = colnames(Y.list[[jj]])
    rownames(ICAp.res.list[[jj]]$Z) = rownames(Y.list[[jj]])
    colnames(ICAp.res.list[[jj]]$Z) = paste0("LV",1:k)
  }#for jj
  
  return(list(ICAp.res.list=ICAp.res.list, B=newB, k = k, number.of.iter = ii, BdiffTrace = BdiffTrace))
}#multi_mirror




simpleDecomp=function(Y, k,svdres=NULL, L1=NULL, L2=NULL,
                      max.iter=200, tol=5e-6, trace=F,
                      rseed=NULL, B=NULL, scale=1,  adaptive.frac=0.05, adaptive.iter=30){
  
  pos.adj=3
  
  ng=nrow(Y)
  ns=ncol(Y)
  
  Bdiff=Inf
  BdiffTrace=double()
  BdiffCount=0
  message("****")
  
  if(is.null(svdres)){
    
    message("Computing SVD")
    set.seed(123);svdres=rsvd(Y, k = k) 
    svdres=rotateSVD(svdres)
  }#if
  
  if(is.null(L1)){
    L1=svdres$d[k]*scale
    if(!is.null(pos.adj)){
      L1=L1/pos.adj
    }
  }
  
  if(is.null(L2)){
    L2=svdres$d[k]*scale
  }
  #    L1=svdres$d[k]/2*scale
  print(paste0("L1 is set to ",L1))
  print(paste0("L2 is set to ",L2))
  
  if(is.null(B)){
    #initialize B with svd
    message("Init")
    B=t(svdres$v[1:ncol(Y), 1:k]%*%diag(sqrt(svdres$d[1:k])))
    #   B=t(svdres$v[1:ncol(Y), 1:k]%*%diag((svdres$d[1:k])))
    #   B=t(svdres$v[1:ncol(Y), 1:k])
  }else{
    message("B given")
  }#else
  
  
  if (!is.null(rseed)) {
    message("using random start")
    set.seed(rseed)
    B = t(apply(B, 1, sample))
  }#if
  
  round2=function(x){signif(x,4)}
  
  getT=function(x){-quantile(x[x<0], adaptive.frac)}
  
  
  for (i in 1:max.iter){
    #main loop    
    Zraw=Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
    
    if(i>=adaptive.iter && adaptive.frac>0){
      
      cutoffs=apply(Zraw,2, getT)
      
      for(j in 1:ncol(Z)){
        Z[Z[,j]<cutoffs[j],j]=0
      }#for j
    }else{
      Z[Z<0]=0
    }#else
    
    oldB=B
    
    B=solve(crossprod(Z)+L2*diag(k))%*%(t(Z)%*%Y)
    
    #update error
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
    err0=sum((Y-Z%*%B)^2)+sum((Z)^2)*L1+sum(B^2)*L2
    if(trace){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }
    
    #check for convergence
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
    }else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }#else if
    
    if(Bdiff<tol &&i>40){
      message(paste0("converged at  iteration ", i))
      break
    }
    if( BdiffCount>5&&i>40){
      message(paste0("stopped at  iteration ", i, " Bdiff is not decreasing"))
      break
    }#if
    
  }#for i
  
  rownames(B)=colnames(Z)=paste("LV",1:k)
  Zproject=Z%*%solve(crossprod(Z)+L2*diag(k))
  return(list(B=B, Z=Z, Zraw=Zraw, Zproject=Zproject,L1=L1, L2=L2))
}#simpleDecomp




#not sure what this is used for
sc_mirror <- function(Y.list, Lap, ICAp.res.list, right.shur.list, L3.list, k, max.iter=200, tol=5e-6, trace = F){
  #Y, Lap and ICAp.res.list are lists
  k <- nrow(ICAp.res.list[[1]]$B)
  number.of.slice <- length(Y.list)
  
  
  #align the genes
  common.genes <- gene.align(Y.list)
  for (ii in 1:number.of.slice){
    Y.list[ii] <- Y.list[[ii]][common.genes,]
    ICAp.res.list[[ii]]$Z <- ICAp.res.list[[ii]]$Z[common.genes, ]
  }#for ii
  
  
  
  BdiffTrace <- c()
  BdiffCount <- 0
  
  
  for (ii in 1:max.iter){
    
    oldB <- lapply(ICAp.res.list, function(x){x$B})
    Z.total <- do.call("+", lapply(ICAp.res.list, function(x){x$Z}))
    
    ######################################################################################################################
    for (jj in 1:length(Y.list)){
      Y <- Y.list[[jj]]
      Z <- ICAp.res.list[[jj]]$Z
      B <- ICAp.res.list[[jj]]$B
      Z0 <- Z.total-ICAp.res.list[[jj]]$Z
      L <- Lap[[jj]]
      L1 <- ICAp.res.list[[jj]]$L1
      L2 <- ICAp.res.list[[jj]]$L2
      L3 <- L3.list[[jj]]
      L4 <- ICAp.res.list[[jj]]$L4
      right.shur <- right.shur.list[[jj]] #already take L4 into consideration
      attr <- ICAp.res.list[[jj]]$attr
      
      #Z update
      ##################
      Zraw=Z=(eigenMapMatMult(Y, t(B))+L3*Z0)%*%solve(tcrossprod(B)+(L1+L3*(number.of.slice-1))*diag(k))
      Z[Z<0]=0
      
      
      #B update
      ##################
      if (attr == "spatial"){
        left <- 2*crossprod(Z)+2*L2*diag(k)
        total <- 2*t(Z) %*% Y
        B <- sylvester_pre(left, right.shur$U, right.shur$S, total)
      }else{
        #single cell case
        B <- solve(crossprod(Z)+L2*diag(k))%*%(t(Z)%*%Y)
      }#else
      
      
      #assign the value back
      ##################
      ICAp.res.list[[jj]]$Z <- Z
      ICAp.res.list[[jj]]$Zraw <- Zraw
      ICAp.res.list[[jj]]$B <- B
      
      #update Z.total
      ##################
      Z.total <- Z0+Z
      
    }#for jj
    ######################################################################################################################
    
    
    #update error
    ####################################################
    #calculate the difference
    Bdiff <- c()
    newB <- lapply(ICAp.res.list, function(x){x$B})
    
    for (jj in 1:length(newB)){
      Bdiff <- c(Bdiff, sum((newB[[jj]]-oldB[[jj]])^2)/sum(oldB[[jj]]^2))
    }#for jj
    
    Bdiff.mean <- mean(Bdiff)
    BdiffTrace <- c(BdiffTrace, Bdiff.mean)
    
    #calculate the error
    ########################
    if (trace & (ii%%1==0) ){
      print(Bdiff.mean)
      #print(BdiffTrace)
    }#if 
    
    #check the convergence
    ########################
    if(ii>52){
      if (Bdiff.mean > BdiffTrace[ii-50]){
        BdiffCount <- BdiffCount+1  
      }else if(BdiffCount>1){
        BdiffCount <- BdiffCount-1
      }#else if 
    }else if(BdiffCount>1){
      BdiffCount <- BdiffCount-1
    }#else if 
    
    #if(Bdiff.mean < tol & ii>40){
    if(Bdiff.mean < tol & ii>20){
      message(paste0("converged at  iteration ", ii))
      break
    }#if
    
    #if( BdiffCount>5 & ii>40){
    if( BdiffCount>5 & ii>20){
      message(paste0("stopped at  iteration ", ii, " Bdiff.mean is not decreasing"))
      break
    }#if
    
  }#for ii
  
  #assign names
  for (jj in 1:length(Y.list)){
    rownames(newB[[jj]]) = rownames(ICAp.res.list[[jj]]$B) = paste0("LV",1:k)
    colnames(newB[[jj]]) = colnames(ICAp.res.list[[jj]]$B) = colnames(Y.list[[jj]])
    rownames(ICAp.res.list[[jj]]$Z) = rownames(Y.list[[jj]])
    colnames(ICAp.res.list[[jj]]$Z) = paste0("LV",1:k)
  }#for jj
  
  return(list(ICAp.res.list=ICAp.res.list, B=newB, k = k, number.of.iter = ii, BdiffTrace = BdiffTrace))
}#sc_mirror





custom.Jaccard <- function(x1, x2){
  sum(x1&x2)/sum(x1|x2)
}#custom.Jaccard


Z.align.jaccard <- function(Z1, Z2, plot = T, num.of.genes = 50){
  align.res <- cor.align.jaccard(Z1, Z2, num.of.genes = num.of.genes)
  
  if (plot){
    cor.res <- align.res$cor.res
    pheatmap::pheatmap(cor.res[,align.res$adjust.index], display_numbers = T, fontsize = 20, cluster_rows = F, cluster_cols = F)
  }else{
    return(align.res$adjust.index)
  }#else
}#Z.align.jaccard


cor.align.jaccard <- function(mat1, mat2, num.of.genes = 50){
  gene.intersect <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[gene.intersect,]
  mat2 <- mat2[gene.intersect,]
  
  #summary(proxy::pr_DB)
  #https://cran.r-project.org/web/packages/proxy/vignettes/overview.pdf
  mat1.top <- apply(mat1,2,function(x){x>=sort(x, decreasing = T)[num.of.genes]})
  mat2.top <- apply(mat2,2,function(x){x>=sort(x, decreasing = T)[num.of.genes]})
  cor.res <- proxy::dist(mat1.top, mat2.top, method = custom.Jaccard, by_rows = F)
  cor.res.output <- cor.res
  
  #fix 1 and find order in 2
  mat1.index <- c()
  mat2.index <- c()
  
  mat1.index.available <- 1:ncol(mat1)
  mat2.index.available <- 1:ncol(mat2)
  
  while (length(mat1.index.available) >0){
    temp <- which(cor.res == max(cor.res), arr.ind = T)
    mat1.index <- c(mat1.index, temp[1,1])
    mat2.index <- c(mat2.index, temp[1,2])
    
    #adjust cor.res and other two variables
    mat1.index.available <- setdiff(mat1.index.available, mat1.index)
    mat2.index.available <- setdiff(mat2.index.available, mat2.index)
    cor.res[temp[1,1],] <- -1
    cor.res[,temp[1,2]] <- -1
  }#while
  
  return(list(cor.res = cor.res.output, adjust.index = mat2.index[order(mat1.index)]))
}#cor.align.jaccard

