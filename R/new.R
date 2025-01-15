#' Adaptive Leiden/Louvain clustering
#'
#' When the number of clusters is set, run louvain/leiden clustering in an adaptive manner.
#' 
#' @export
leiden_adaptive <- function(nn, num.of.cluster = 5, resolution.start = 0.5, adaptive.size = 2, method = "louvain", verbose = T){
  
  #partition <- leiden::leiden(nn, resolution_parameter = resolution.start)
  if (method == "leiden"){
    partition <- Seurat:::RunLeiden(nn, resolution.parameter = resolution.start, method = "matrix") 
  }else if (method =="louvain"){
    partition <- Seurat:::RunModularityClustering(nn, resolution = resolution.start, algorithm = 1, print.output = F) 
  }#else if
  
  
  if (length(unique(partition)) < num.of.cluster){
    reso.left <- resolution.start
    reso.right <- resolution.start*adaptive.size
    reso <- reso.right
  }else if (length(unique(partition)) > num.of.cluster){
    reso.left <- resolution.start/adaptive.size
    reso.right <- resolution.start
    reso <- reso.left
  }else{
    return(partition)
  }#else
  
  
    
    while(length(unique(partition)) != num.of.cluster){
 
      if (reso < 0){
        stop("The resolution parameter is not valid.")
      }#if
      
      if (abs(reso.left-reso.right) < 1e-4){
        message("Re-initalize.") 
        return( leiden_adaptive(nn, num.of.cluster, resolution.start = resolution.start+0.1, adaptive.size, method, verbose) )
        
      }else{
      
        if (method == "leiden"){
          #partition <- leiden::leiden(nn, resolution_parameter = reso)
          partition <- Seurat:::RunLeiden(nn, resolution.parameter = reso, method = "matrix")  
        }else if (method == "louvain"){
          partition <- Seurat:::RunModularityClustering(nn, resolution = reso, algorithm = 1, print.output = F)  
        }#else if
        
        if ( length(unique(partition)) < num.of.cluster){
          message("Explode")
          if (reso == reso.left){
            reso.left <- reso.left
            reso <- (reso.left+reso.right)/2
            reso.right <- reso.right
            
          }else if (reso < reso.right){
            reso.left <- reso
            reso <- (reso+reso.right)/2
            reso.right <- reso.right
          }else if (reso == reso.right){
            reso <- reso.right*adaptive.size
            reso.left <- reso.right
            reso.right <- reso
          }#else if
          
        }else if ( length(unique(partition)) > num.of.cluster){
          message("Shrink")
          
          if (reso == reso.left){
            reso <- reso/adaptive.size
            reso.right <- reso.left
            reso.left <- reso
            
          }else if (reso < reso.right){
            reso.right <- reso
            reso <- (reso.left+reso)/2
            reso.left <- reso.left
            
          }else if (reso == reso.right){
            reso.left <- reso.left
            reso <- (reso.left+reso.right)/2
            reso.right <- reso.right
            
          }#else if
          
        }else{
          return(partition)
        }#else
        
        if (verbose){
          cat(paste0(reso, " ", reso.left, " ", reso.right, "\n"))  
        }#if verbose  
        
      }#else
      
    }#while
    
 
}#leiden_adaptive






  
  self_extract <- function(proj, ids, opt.in = "Z"){
    
    res <- c()
    
    ii.list <- as.numeric(unlist(lapply(strsplit(ids, "_"), function(x){x[1]})))
    jj.list <- as.numeric(unlist(lapply(strsplit(ids, " "), function(x){x[2]})))
    
    if (opt.in=="Z"){
      #extract B      
      for (ii in 1:length(ii.list)){
        res <- rbind(res, proj[[ii.list[ii]]]$B[jj.list[ii],])
      }#for ii
      rownames(res) <- ids
      colnames(res) <- colnames(proj[[1]]$B)
      
    }else if (opt.in == "B"){
      #extract Z
      for (ii in 1:length(ii.list)){
        res <- rbind(res, proj[[ii.list[ii]]]$Z[,jj.list[ii]])
      }#for ii
      rownames(res) <- ids
      colnames(res) <- rownames(proj[[1]]$Z)

    }#else if
    
    return(res)
    
  }#self_extract
  
  

#' Consensus Decomposition
#'
#' Figure out the consensus of latent variables
#' 
#' @export
self_deco <- function(proj,  LVs.filter.thr = 0.9, freq = 1, opt = "B"){
  
  num.of.slice <- length(proj)
  
  #row names and colnames
  row.names <- c()
  LVs <- c()
  
  if (opt == "B"){
    
    col.names <- colnames(proj[[1]]$B)
    
    for (ii in 1:num.of.slice){
      LVs <- rbind(LVs, proj[[ii]]$B)
      row.names <- c(row.names, paste0(ii, "_", rownames(proj[[ii]]$B)) )
      if (!all(col.names==colnames(proj[[ii]]$B))){
        warning("column samples are mismatched at ", ii)
        break
      }#if
    }#for ii
    

  }else if (opt == "Z"){
    col.names <- rownames(proj[[1]]$Z)
    
    for (ii in 1:num.of.slice){
      LVs <- rbind(LVs, t(proj[[ii]]$Z))
      row.names <- c(row.names, paste0(ii, "_", colnames(proj[[ii]]$Z)) )
      if (!all(col.names==rownames(proj[[ii]]$Z))){
        warning("column samples are mismatched at ", ii)
        break
      }#if
    }#for ii

  }#else if opt == "Z"
  
  
  #construct the LVs
  ###########################################
  
  rownames(LVs) <- row.names
  colnames(LVs) <- col.names
    
  
  #filter by correlation
  ###########################################
  
       
      #only selected the most highly correlated ones
      cor.res <- cor(t(LVs))
      
      #converted to a connected graph based on the threshold
      adj_matrix <- cor.res > LVs.filter.thr
      adj_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = F)
      components_list <- igraph::components(adj_graph)
      components_selected <- which(components_list$csize>freq)
      components_selected <- components_selected[order(components_list$csize[components_selected] , decreasing = T)]
      
      print(sort(components_list$csize, decreasing = T))
      
      #output the module size
      components.size <- components_list$csize[components_selected]
          
      #only pick up one from correlated blocks
      include.names <- c()
      for (ii in 1:length(components_selected)){
        #select the one with smaller index
        include.names <- c(include.names, names(components_list$membership)[which(components_list$membership==components_selected[ii])][1] )
      }#for ii
      
      LVs <- LVs[include.names,]



#extract the corresponding ones
  if (opt == "B"){
    LVs.pair <- self_extract(proj, rownames(LVs), opt.in = "B")
  }else if (opt == "Z"){
    LVs.pair <- self_extract(proj, rownames(LVs), opt.in = "Z")
  }#else if
  
  return(list(LVs=LVs, LVs.pair = LVs.pair, components.size = components.size))
}#self_deco




