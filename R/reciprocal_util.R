

GenePlot <- function(dat, features, pt.size = 2, ratio = 1, use.myratio = F){
  coord <- Seurat::GetTissueCoordinates(object = dat@images[[1]])
  myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  message(paste0("The aspect ratio is ", myratio))
  
  if (use.myratio){
    ratio <- myratio
  }#if
  
  SpatialFeaturePlot(dat, features = features, pt.size.factor = pt.size)+ggplot2::theme(aspect.ratio = ratio) 
}#GenePlot




#input should be the gene, L and ICAp.res
reciprocal_default <- function(Y.list, Lap, ICAp.res.list, num.of.reference.slice = NULL){
  
  num.of.slice <- length(Y.list)
  if (is.null(num.of.reference.slice)){
    num.of.reference.slice <- length(ICAp.res.list)  
  }#if
  #else it will be first several slices based on the provided numbers
  message(paste0("The number of refernce slices is: " , num.of.reference.slice))
  message(paste0("The number of total slices is: ", num.of.slice))

  #calculate the projection and construct a 2d list
  ##############################################################################
  proj <- list()
  
  for (ii in 1:num.of.reference.slice){
    proj.ii.on.jj <- list()
    
    for (jj in 1:num.of.slice){
    
      if (ii == jj){
        proj.ii.on.jj[[jj]] <- ICAp.res.list[[jj]]
        
      }else{
          proj.ii.on.jj[[jj]] <- impact(ICAp.res.list[[ii]]$Z, Y.list[[jj]], Lap[[jj]], ICAp.res.list[[jj]]$L4, query.shur0 = ICAp.res.list[[jj]]$shur0, query.ICAp.res = ICAp.res.list[[jj]], scale= 1, max.iter = 200, cor.thr = 0.8)
          
      }#else
      
    }#for jj 
    
    proj[[ii]] <- proj.ii.on.jj
  }#for ii
  
  return(proj)
  
}#reciprocal_default


p.intersect <- function(x){
  
  res <- x[[1]]
  if (length(x)>1){
    
    for (ii in 2:length(x)){
      res <- intersect(res, x[[ii]])
    }#for ii
  }#if
  
  return(res)
}#p.intersect



#' apply the existing integration result on the new dataset
#' @export
reciprocal_with_Z <- function(Y.list, Lap, ICAp.res.list, ICAp.res.list.reference, Z.names){
  
  num.of.slice <- length(Y.list)
  num.of.reference.slice <- length(ICAp.res.list.reference)
  message(paste0("The number of refernce slices is: " , num.of.reference.slice))
  message(paste0("The number of total slices is: ", num.of.slice))
  
  
  #first construct Z
  ##########################################################
  slice.index <- as.numeric(lapply(strsplit(Z.names, "_"), function(x){x[1]}))
  LV.names <- as.character(lapply(strsplit(Z.names, "_"), function(x){x[2]}))
  
  Z.inuse.list <- list()
  for (ii in 1:length(Z.names)){
      Z.inuse.list[[ii]] <- ICAp.res.list.reference[[slice.index[ii]]]$Z[,LV.names[ii]]
  }#for ii
  
  gene.names <- p.intersect(lapply(Z.inuse.list, names))
    
  Z.inuse <- do.call(cbind, lapply(Z.inuse.list, function(x){x[gene.names]}))
  rownames(Z.inuse) <- gene.names
  colnames(Z.inuse) <- Z.names
  
  #calculate the projection
  ##########################################################
  proj <- list()
  
  for (ii in 1:num.of.slice){
    proj[[ii]] <- impact(Z.inuse, Y.list[[ii]], Lap[[ii]], ICAp.res.list[[ii]]$L4, query.shur0 = ICAp.res.list[[ii]]$shur0, query.ICAp.res = ICAp.res.list[[ii]], scale= 1, max.iter = 200, cor.thr = 0.8)
  }#for ii
  
  LVs.inuse <- do.call(cbind, lapply(proj, function(x){x$B}))
  
  col.names <- c()
  for (ii in 1:num.of.slice){
    col.names <- c(col.names, paste0(ii, "_", colnames(Y.list[[ii]])) )
  }#for ii
  
  rownames(LVs.inuse) <- Z.names
  colnames(LVs.inuse) <- col.names
  
  return(LVs.inuse)
}#reciprocal_with_Z







#construct the B.list and global LVs with different options
#' Finalize the integration by filtering out duplicated LVs based on similarity and strategy
#' @export
reciprocal_deco <- function(proj, option = "L2norm", LVs.filter = F, LVs.filter.thr = 0.9, mod = "all"){
  
  num.of.slice <- length(proj)
  
  #row names and colnames
  row.names <- c()
  col.names <- c()
  for (ii in 1:num.of.slice){
    row.names <- c(row.names, paste0(ii, "_", rownames(proj[[ii]][[ii]]$B)) )
    col.names <- c(col.names, paste0(ii, "_", colnames(proj[[ii]][[ii]]$B)) )
  }#for ii
  
  
  #merge
  dim.list <- list()
  
  for (ii in 1:num.of.slice){
    if (option == "default"){
      dim.list[[ii]] <- do.call(cbind, lapply(proj[[ii]], function(x){x$B}) )
    }else if (option == "L2norm"){
      dim.list[[ii]] <- do.call(cbind, lapply( lapply(proj[[ii]], function(x){x$B}), L2Norm))
    }else if (option == "L2norm.joint"){
      dim.list[[ii]] <- L2Norm(do.call(cbind, lapply(proj[[ii]], function(x){x$B}) ))
    }#else if
  }#for ii
  
  #construct the LVs
  ###########################################
  LVs <- do.call(rbind, dim.list)
  
  rownames(LVs) <- row.names
  colnames(LVs) <- col.names
  
  
  
  #filter by correlation
  ###########################################
  
  if (LVs.filter){
  
    
    if (mod == "all"){
      
      #use as many LVs < LVs.filter.thr as possible to generate the final LVs in use
      
      cor.res <- cor(t(LVs))
      diag(cor.res) <- 0
      cor.res[lower.tri((cor.res))] <- 0
      #iteratively remove the correlations
      while (max(cor.res)> LVs.filter.thr){
        #remove the larger LV index
        LV.to_drop <- max(which(cor.res==max(cor.res[upper.tri(cor.res)]), arr.ind = T))
        LVs <- LVs[-LV.to_drop,]
        
        cor.res <- cor(t(LVs))
        diag(cor.res) <- 0
        cor.res[lower.tri((cor.res))] <- 0
      }#while 
   
         
    }else if (mod == "common"){
      
      #only selected the most highly correlated ones
      cor.res <- cor(t(LVs))
      
      #converted to a connected graph based on the threshold
      adj_matrix <- cor.res > LVs.filter.thr
      adj_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = F)
      components_list <- igraph::components(adj_graph)
      components_selected <- which(components_list$csize>1)
      
      #only pick up one from correlated blocks
      include.names <- c()
      for (ii in 1:length(components_selected)){
        #select the one with smaller index
        include.names <- c(include.names, names(components_list$membership)[which(components_list$membership==components_selected[ii])][1] )
      }#for ii
      
      LVs <- LVs[include.names,]
      
    }#else if mod == "common"
    
  }#if LVs.filter
  
  return(LVs)
}#reciprocal_deco











