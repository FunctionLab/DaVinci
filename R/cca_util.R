
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
#' For scalable implementation
#' @export
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
#' Used by reciprocal_default() and reciprocal_with_Z()
impact <- function(reference.Z, query.gene.exp, query.L, query.L4, query.shur0 = NULL, query.ICAp.res = NULL, scale=1, max.iter = 200, cor.thr = 0.8, save.complete=F){
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
    return(list(B=B, L1 = L1, L2 = L2, L4 = L4, k = k))
    #return(list(B=B, L1 = L1, L2 = L2, L4 = L4, k = k, query.shur0 = query.shur0, right.shur = right.shur))
  }#else
  
}#impact







#' prepare the input data
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





#' prepare the input data
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









