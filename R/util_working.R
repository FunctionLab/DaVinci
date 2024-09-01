require(cluster) #silhouette function
#Sys.getenv("LD_LIBRARY_PATH"), path is not included
#dyn.load("/mnt/home/wmao1/software/miniconda3/lib/libgdal.so.20")
#library(terra)
require(emojifont) #recommended but not required






#' Automatically tune the number of latent variables
#' @importFrom rsvd rsvd
#' @export
AutoTune <- function(gene.exp, L, k.arg.list = c(5, 10, 12, 15, 18, 20, 25), L4.arg.list=exp(seq(log(100), log(10), length.out = length(k.arg.list))), frac.impute=0.1, frac.seed = 1 ){
  #create the dataset
  nonzero.index <- which(gene.exp!=0)
  set.seed(frac.seed)
  sample.index <- sample(nonzero.index, length(nonzero.index)*frac.impute)
  gene.exp.hide <- gene.exp
  gene.exp.hide[sample.index] <- 0

  recon.grid <- c()
  L4.grid <- c()


  #calculate SVD ahead of time to save some time
  #########################
  message("Computing SVD")
  set.seed(123)
  svdres=rsvd::rsvd(gene.exp.hide, k = max(k.arg.list)+1) #to_drop -> +1
  svdres=rotateSVD(svdres)
  #########################


  for (ii in 1:length(k.arg.list)){
    print(ii)

    ICAp.res.impute <- manifoldDecomp_adaptive(gene.exp.hide, L, k = k.arg.list[ii], L4 = L4.arg.list[ii], L4_adaptive = 2, to_drop = T, svdres = svdres)

    #reconstruction error
    Y.recon <- ICAp.res.impute$Z %*% ICAp.res.impute$B
    recon.grid <- c(recon.grid, recon.error(gene.exp, Y.recon, sample.index))
    L4.grid <- c(L4.grid, ICAp.res.impute$L4)
  }#for ii



  #scan for local optimal
  #if there are many, pick up the larger one
  ##########################################################
  temp <- c()
  temp.seq <- 2: (length(recon.grid)-1)
  for (kk in  temp.seq){
    if ( ( recon.grid[kk] <=recon.grid[kk-1]) & (recon.grid[kk] <=recon.grid[kk+1]) ){
      temp <- c(temp, kk)
    }#if
  }#for ii
  k.to_report <- k.arg.list[max(temp)]


  #return the result
  ################################
  return(list(k.arg.list = k.arg.list, L4.arg.list = L4.arg.list, recon.grid = recon.grid, L4.grid = L4.grid, k.optim = k.to_report))
}#AutoTune



#' Return the ARI and NMI values of clusters
#' @export
ClusterMetric <- function(x, y){
    include.index <- which(!is.na(x) & !is.na(y))
    message(paste0("ARI: ", aricode::ARI(x[include.index], y[include.index])))
    message(paste0("NMI: ", aricode::NMI(x[include.index], y[include.index])))
}#ClusterMetric




unicode.convert <- function(x){
  switch(x,
         "1" = "\u2460",
         "2" = "\u2461",
         "3" = "\u2462",
         "4" = "\u2463",
         "5" = "\u2464",
         "6" = "\u2465",
         "7" = "\u2466",
         "8" = "\u2467",
         "9" = "\u2468",
         "10" = "\u2469",
         "11" = "\u246A",
         "12" = "\u246B",
         "13" = "\u246C",
         "14" = "\u246D",
         "15" = "\u246E",
         "16" = "\u246F",
         "17" = "\u2470",
         "18" = "\u2471",
         "19" = "\u2472",
         "20" = "\u2473")
}#unicode.convert

#scCustomize



#https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' Return the color palette given cluster vector
#' @export
colorPalette <- function(partition, fill = T){
  num.of.clusters <- length(unique(partition))
  val.names <- sort(as.character(1:num.of.clusters))

  if (num.of.clusters <=8){
    #ColorBlind_Pal()
    values <- ggthemes::colorblind_pal()(num.of.clusters)
  }else if (num.of.clusters <=36){
    values <- Polychrome::palette36.colors()[1:num.of.clusters]
  }else{
    values <- scales::hue_pal()(num.of.clusters)
  }#else

  #names(values) <- val.names
  names(values) <- as.character(unique(partition))
  if (fill){
    return(ggplot2::scale_fill_manual(values = values))
  }else{
    return(ggplot2::scale_color_manual(values = values))
  }#else

}#colorPalette





#caveat: cluster may be not the real cluster: Seurat treat cluster as
#' Plot clusters on top of Seurat object
#' For Seurat object
#' @details This function requires `Seurat`. Make sure it is installed
#' @export
ClusterPlot <- function(dat, cluster.label, pt.size = 2, ratio = 1, use.myratio = F, cluster.highlight = NULL, image.alpha = 1){
  dat.tmp <- dat

if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for this function but is not installed. Please install it.")
  }
  #coord <- Seurat::GetTissueCoordinates(object = dat@images$slice1)

  if (use.myratio){
    coord <- Seurat::GetTissueCoordinates(object = dat@images[[1]])
    myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
    message(paste0("The aspect ratio is ", myratio))

    ratio <- myratio
  }#if

  dat.tmp@meta.data$Cluster <- cluster.label
  if (is.null(cluster.highlight)){
    SpatialPlot(dat.tmp, group.by = "Cluster", pt.size.factor = pt.size, image.alpha = image.alpha)+ggplot2::theme(aspect.ratio = ratio)
  }else{
    SpatialPlot(dat.tmp,  pt.size.factor = pt.size, cells.highlight = rownames(dat.tmp@meta.data)[dat.tmp$Cluster==cluster.highlight], cols.highlight = c("red",ggplot2::alpha("gray10", 0)),alpha = NULL, image.alpha = image.alpha)+ggplot2::theme(aspect.ratio = ratio)
  }#else

}#ClusterPlot




#' Plot clusters on top of Seurat object, display cluster lables per spot
#' For Seurat object
#' @import ggpubr
#' @details This function requires `Seurat`. Make sure it is installed
#' @export
ClusterPlot_label <- function(dat, cluster.label, pt.size = 2, text.size = 7, ratio = 1, use.myratio = F, image.alpha = 1){
  dat.tmp <- dat

if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for this function but is not installed. Please install it.")
  }

  if (use.myratio){
    coord <- Seurat::GetTissueCoordinates(object = dat@images[[1]])
    myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
    message(paste0("The aspect ratio is ", myratio))

    ratio <- myratio
  }#if

  dat.tmp@meta.data$Cluster <- cluster.label
  p1 <- SpatialPlot(dat.tmp, group.by = "Cluster", pt.size.factor = pt.size, image.alpha=image.alpha)+ggplot2::theme(aspect.ratio = ratio)

  data.p1 <- ggplot_build(p1)
  #data.p1$data[[1]]$shape <- NA #not necessary
  dat.extra <- data.p1$data[[1]]

  grp.real <- sort(unique(cluster.label))
  names(grp.real) <- sort(unique(as.character(dat.extra$group)))


  dat.extra$label <- unlist(lapply(grp.real[dat.extra$group], unicode.convert))
  dat.extra$Cluster <- grp.real[dat.extra$group]

  p1+geom_text(data = dat.extra, aes(x=x,y=y, label = label), size = text.size)
}#ClusterPlot_label



#' Plot clusters on top of Seurat object, no background image
#' For Seurat object
#' @import ggpubr
#' @details This function requires `Seurat`. Make sure it is installed
#' @export
ClusterPlot_fast <- function(dat, cluster.label, pt.size=2, text.size=7, ratio = 1, use.myratio=F, label.on = F, cluster.highlight = NULL){

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for this function but is not installed. Please install it.")
  }

  if (use.myratio){
    coord <- Seurat::GetTissueCoordinates(object = dat@images[[1]])
    myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
    message(paste0("The aspect ratio is ", myratio))

    ratio <- myratio
  }#if

  #dat.test <- data.frame(x=dat@images$slice1@coordinates$row, y = dat@images$slice1@coordinates$col, cluster = cluster.label)
  #flip the coordinate
  dat.test <- data.frame(x=dat@images[[1]]@coordinates$imagecol, y = dat@images[[1]]@coordinates$imagerow, cluster = cluster.label)

  if (is.null(cluster.highlight)){

    if (label.on){
      dat.test$label <- unlist(lapply(dat.test$cluster, unicode.convert))
      ggplot(dat.test, aes(x=x,y=-y,color= cluster))+geom_text(aes(label=label), size= text.size)+theme_void()+ggplot2::theme(aspect.ratio = ratio)+colorPalette(dat.test$cluster, fill = F)

    }else{
      ggplot(dat.test, aes(x=x,y=-y,color= cluster))+geom_point(size= pt.size)+theme_void()+ggplot2::theme(aspect.ratio = ratio)+colorPalette(dat.test$cluster, fill = F)
    }#else

  }else{
    #highlight one single cluster
    cluster.unique <- unique(dat.test$cluster)
    values <- rep("grey", length(cluster.unique))
    names(values) <- cluster.unique
    values[as.character(cluster.highlight)] <- "red"

    ggplot(dat.test, aes(x=x,y=-y,color= cluster))+geom_point(size= pt.size)+theme_void()+ggplot2::theme(aspect.ratio = ratio)+ggplot2::scale_color_manual(values = values)
  }#else

}#ClusterPlot_fast






#' Visualize discrete variables only with coordinate information
#' @import ggpubr
#' @export
scatter.DiscretePlot <- function(coor, cluster.label, pt.size = 2, ratio = NULL, plot.all = F){

  cluster.label <- as.character(cluster.label)
  cluster.unique <- sort(unique(cluster.label))

  dat.plot <- data.frame(x= coor[,1], y = coor[,2], cluster = cluster.label)

   if (plot.all){

    p <- list()
    for (ii in 1:length(cluster.unique)){
      cluster.highlight <- cluster.unique[ii]

      values <- rep("grey", length(cluster.unique))
      names(values) <- cluster.unique
      values[as.character(cluster.highlight)] <- "red"

      p[[ii]] <- ggplot(dat.plot, aes(x=x,y=y,color= cluster))+geom_point(size= pt.size)+theme_void()+ggplot2::scale_color_manual(values = values)
      if (!is.null(ratio)){
          p[[ii]] <- p[[ii]]+ggplot2::theme(aspect.ratio = ratio)
      }#if

    }#for ii

    p <- ggarrange(plotlist = p)

  }else{

    p <- ggplot(dat.plot, aes(x=x,y = y, col = cluster))+geom_point(size = pt.size)+theme_void()+colorPalette(dat.plot$cluster, fill = F)
    if (!is.null(ratio)){
          p <- p+ggplot2::theme(aspect.ratio = ratio)
    }#if

  }#else

  return(p)
}#scatter.DiscretePlot 

#DiscretePlot <- function(coor, cluster.label, pt.size=2){
  #dat.plot <- data.frame(x=coor[,1], y = coor[,2], cluster = cluster.label)
  #ggplot(dat.plot, aes(x=x,y=-y,color= cluster))+geom_point(size= pt.size)+theme_void()+colorPalette(dat.plot$cluster, fill = F)
#}#DiscretePlot






#seurat object plot
#' Visualize single latent variable
#' For Seurat object
#' @import ggpubr
#' @details This function requires `Seurat`. Make sure it is installed
#' @export
LvPlot <- function(dat, loading, LVs, LV.index = 1, gene.verbose =T, pt.size = 2, ratio = 1, use.myratio = F, ws = F, x.offset = 5, image.alpha= 1){
  dat.tmp <- dat

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for this function but is not installed. Please install it.")
  }

  if (use.myratio){
    #coord <- Seurat::GetTissueCoordinates(object = dat@images$slice1)
    coord <- Seurat::GetTissueCoordinates(object = dat@images[[1]])
    myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
    message(paste0("The aspect ratio is ", myratio))

    ratio <- myratio
  }#if

  print(paste0("LV ", LV.index))

  if (gene.verbose){
    print(sort(loading[,LV.index], decreasing = T)[1:10])
  }#if

  tmp <- LVs[LV.index,]
  dat.tmp@meta.data$LV <- tmp
  if (length(dat.tmp@images[[1]]@key)==0){
    dat.tmp@images[[1]]@key <- "image_"
  }#if

  res <- Seurat::SpatialPlot(dat.tmp, features = "LV", pt.size.factor = pt.size, image.alpha=image.alpha)+ggplot2::theme(aspect.ratio = ratio)

  if (gene.verbose){
    res <- res+annotate("text", x = layer_scales(res)$x$range$range[2]+x.offset, y = median(layer_scales(res)$y$range$range), label = paste(names(sort(loading[,LV.index], decreasing = T)[1:10]), collapse = "\n"), size =3)
  }#gene.verbose

  #add title
  res <- res+ggtitle(paste0("LV ", LV.index))

  #use customer color palette
  res <- res+ggthemes::scale_fill_gradient2_tableau(palette = "Orange-Blue Diverging", trans = "reverse")

  if (ws){
    res <- res+coord_cartesian(clip = "off")
  }#if

  return(res)
}#LvPlot



#' Visualize single latent variable, no background image
#' For Seurat object
#' @import ggpubr
#' @details This function requires `Seurat`. Make sure it is installed
#' @export
LvPlot_fast <- function(dat, val, pt.size=2, ratio = 1, use.myratio=F, plot.all = F){

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for this function but is not installed. Please install it.")
  }

  if (use.myratio){
    coord <- Seurat::GetTissueCoordinates(object = dat@images[[1]])
    myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
    message(paste0("The aspect ratio is ", myratio))

    ratio <- myratio
  }#

  if (!plot.all){
    #flip the coordinate
    dat.test <- data.frame(x=dat@images[[1]]@coordinates$imagecol, y = dat@images[[1]]@coordinates$imagerow, val = val)

    ggplot(dat.test, aes(x=x,y=-y,color= val))+geom_point(size= pt.size)+theme_void()+ggplot2::theme(aspect.ratio = ratio)+ggthemes::scale_color_gradient2_tableau(palette = "Orange-Blue Diverging", trans = "reverse")

  }else{
    message(paste0("Total number of panels is ", nrow(val)))

    plot.list <- list()
    for (ii in 1:nrow(val)){
      dat.test <- data.frame(x=dat@images[[1]]@coordinates$imagecol, y = dat@images[[1]]@coordinates$imagerow, val = val[ii,])
      plot.list[[ii]] <-     ggplot(dat.test, aes(x=x,y=-y,color= val))+geom_point(size= pt.size)+theme_void()+ggplot2::theme(aspect.ratio = ratio)+ggthemes::scale_color_gradient2_tableau(palette = "Orange-Blue Diverging", trans = "reverse")+ggtitle(paste0("LV ", ii))+theme(plot.title = element_text(size= 20, face = "bold"))

    }#for ii
    ggarrange(plotlist = plot.list)
  }#else

}#LvPlot_fast






#' Visualize continous variables only with coordinate information
#' @import ggpubr
#' @export
scatter.FeaturePlot <- function(coor, loading=NULL, LVs, LV.index = 1, gene.verbose = F, pt.size = 2, ratio = NULL, x.offset = 1.05, font.size = 7, plot.all = F){

  if (gene.verbose & is.null(loading)){
    stop("loading can't be NULL if gene names will be printed out")
  }#if

  if (is.numeric(LVs) & is.vector(LVs)){
    #no gene verbose in this case

    dat.plot <- data.frame(x= coor[,1], y = coor[,2], val = LVs)
    p <- ggplot2::ggplot(dat.plot, aes(x=x,y = -y, col = val))+geom_point(size = pt.size)+theme_void()+paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1)
    
    if (!is.null(ratio)){
        p <- p+ggplot2::theme(aspect.ratio = ratio)
    }#if

  }else{
    if (nrow(LVs) > ncol(LVs)){
      LVs <- t(LVs)
    }#if

  #plot all rows
  if (plot.all){
    #no gene.verbose in this case
    p <- list()

    for (ii in 1:nrow(LVs)){
      dat.plot <- data.frame(x= coor[,1], y = coor[,2], val = LVs[ii,])
      p[[ii]] <- ggplot2::ggplot(dat.plot, aes(x=x,y = y, col = val))+geom_point(size = pt.size)+theme_void()+paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1)+ggtitle(paste0("LV ", ii))+ggplot2::theme(aspect.ratio = ratio)
    }#for ii

    p <- ggarrange(plotlist = p)

  }else{
    #plot one specific row
    print(paste0("LV ", LV.index))

    if (gene.verbose){
      print(sort(loading[,LV.index], decreasing = T)[1:10])
    }#if

    dat.plot <- data.frame(x= coor[,1], y = coor[,2], val = LVs[LV.index,])

    p <- ggplot(dat.plot, aes(x=x,y = y, col = val))+geom_point(size = pt.size)+theme_void()+paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1)+ggtitle(paste0("LV ", LV.index))+ggplot2::theme(aspect.ratio = ratio)

    if (gene.verbose){
      p <- p+annotate("text", x = max(dat.plot$x)*x.offset , y = median(dat.plot$y), label = paste(names(sort(loading[,LV.index], decreasing = T)[1:15]), collapse = "\n"), size = font.size)
    } #if
  }#else

  }#else, matrix input
  
  return(p)
}#scatter.FeaturePlot

#FeaturePlot <- function(coor, val, pt.size=2){
  #dat.plot <- data.frame(x=coor[,1], y = coor[,2], val = val)
  #ggplot(dat.plot, aes(x=x,y=-y,color= val))+geom_point(size= pt.size)+theme_void()+ggthemes::scale_color_gradient2_tableau(palette = "Orange-Blue Diverging", trans = "reverse")
#}#FeaturePlot






silhouette_score <- function(input, k, random.seed = 1){
  set.seed(random.seed)
  km <- kmeans(input, centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(input))
  mean(ss[, 3])
}#silhouette_score




leiden_adaptive <- function(nn, num.of.cluster = 5, resolution.start = 0.5, adaptive.size = 2){
  partition <- leiden::leiden(nn, resolution_parameter = resolution.start)

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

  cat(paste0(reso, " ", reso.left, " ", reso.right, "\n"))

  while(length(unique(partition)) != num.of.cluster){
    partition <- leiden::leiden(nn, resolution_parameter = reso)

    if (length(unique(partition)) < num.of.cluster){
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

    cat(paste0(reso, " ", reso.left, " ", reso.right, "\n"))
  }#while

}#leiden_adaptive










md.default <- function(Y, B, k, max.iter, L1, L2, L4, right.shur, adaptive.iter, adaptive.frac, trace, tol, pos, thr = 0.8, early_flag = T){

  Bdiff=Inf
  BdiffTrace=double()
  BdiffCount=0

  getT=function(x){-quantile(x[x<0], adaptive.frac)}
  round2=function(x){signif(x,4)}


  for ( i in 1:max.iter){
    #Z update
    ######################################
    Zraw=Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))

    if (pos){
      if(i>=adaptive.iter && adaptive.frac>0){
        cutoffs=apply(Zraw,2, getT)
        for(j in 1:ncol(Z)){
          Z[Z[,j]<cutoffs[j],j]=0
        }#for j
      }else{
        Z[Z<0]=0
      }#else
    }#if pos



    #B update
    ######################################
    oldB=B
    #B=solve(crossprod(Z)+L2*diag(k))%*%(t(Z)%*%Y)
    left <- 2*crossprod(Z)+2*L2*diag(k)
    total <- 2*t(Z) %*% Y

    B <- sylvester_pre(left, right.shur$U, right.shur$S, total)


    #update error
    ######################################
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
    }#if

    if( BdiffCount>5&&i>40){
      message(paste0("stopped at  iteration ", i, " Bdiff is not decreasing"))
      break
    }#if


    #stop early
    if (early_flag){
      if (i > 20){
        B.cor.res <- cor(t(B))
        B.cor.val <- B.cor.res[upper.tri(B.cor.res)]

        if (sum(B.cor.val>=thr)>1){
          return(list(flag = "Shrink"))
        }#if
      }#if i>20

    }#if

  }#for i


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


  return(list(flag = flag, B = B, Z =Z, Zraw = Zraw))
}#md.default





#' Main function
#' @export
manifoldDecomp_adaptive=function(Y, L, k, svdres=NULL, L1=NULL, L2=NULL, L4 = NULL, shur0 = NULL, max.iter=200, tol=5e-6, trace=F,rseed=NULL, B=NULL, scale=1,  adaptive.frac=0.05, adaptive.iter=30, L4_adaptive =2, to_drop =T, pos = T, cor.thr = 0.8, save.complete = F, verbose = T){

  pos.adj=3
  ng=nrow(Y)
  ns=ncol(Y)

  if (to_drop){
    k <- k+1

    if ( !is.null(B)){

      if ((nrow(B)<k)){
        if (verbose){
          message("Add random LV")
        }#if vebose

        #set.seed(rseed)
        #B.base <- matrix(rnorm(k*ncol(B)), nrow =k, ncol = ncol(B) )
        B.base <- t(matrix(apply(B,2,mean), ncol = k, nrow = ncol(B)))
        B.base[1:nrow(B),] <- B
        B <- B.base
      }#if nrow(B)<k
    }#!is.null B
  }#if to _drop

  if (verbose){
    message("****")
  }#if verbose


  if(is.null(svdres)){

    if (verbose){
      message("Computing SVD")
    }#if verbose

    set.seed(123)
    svdres=rsvd(Y, k = k)

    svdres=rotateSVD(svdres)

    #  show(svdres$d[k])
  }#svdres


  #L1
  #######################################################
  if(is.null(L1)){
    L1=svdres$d[k]*scale
    if(!is.null(pos.adj)){
      L1=L1/pos.adj
    }#if
  }#if is.null(L1)


  #L2
  #######################################################
  if(is.null(L2)){
    L2=svdres$d[k]*scale
  }#if is.null(L2)

  #L1=svdres$d[k]/2*scale
  if (verbose){
    print(paste0("L1 is set to ",L1))
    print(paste0("L2 is set to ",L2))
  }#if verbose



  #B
  #######################################################
  if(is.null(B)){
    #initialize B with svd
    if (verbose){
      message("Init")
    }#if verbose

    B=t(svdres$v[1:ncol(Y), 1:k]%*%diag(sqrt(svdres$d[1:k])))
  }else{
    if (verbose){
      message("B given")
    }#if verbose

  }#B initialization


  if (!is.null(rseed)) {
    if (verbose){
      message("using random start")
    }#if verbose

    set.seed(rseed)
    B = t(apply(B, 1, sample))
  }#is.null rseed



  B0 <- B
  if (is.null(shur0)){
    shur0 <- rcpp_shur(L)
  }else{
    if (verbose){
      message("shur0 has benn set.")
    }#if verbose

  }#else


  #updates: md.default
  #######################################################

  #right <- L4*t(L)+L4*L
  #right <- 2*L4*L
  #right.shur <- rcpp_shur(right)
  right.shur <- shur0
  right.shur$S <- 2*L4*right.shur$S
  if (verbose){
    message("Shur done")
  }#if verbose


  md.run <- md.default(Y, B0, k, max.iter, L1, L2, L4, right.shur, adaptive.iter, adaptive.frac, trace, tol, pos, thr = cor.thr)


  if (md.run$flag == "Explode"){

    if (verbose){
      message("Explode")
    }#if vebose


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
    right.shur <- shur0
    right.shur$S <- 2*L4_pointer*right.shur$S

    #if L4_left and L4_right is too close, then Done
    if (abs(L4_left-L4_right) < 0.1){
      md.run <- md.default(Y, B0, k, max.iter, L1, L2, L4 = L4_pointer, right.shur, adaptive.iter, adaptive.frac, trace, tol, pos, thr = cor.thr, early_flag = F)
      md.run$flag <- "Done"
    }else{
      md.run <- md.default(Y, B0, k, max.iter, L1, L2, L4 = L4_pointer, right.shur, adaptive.iter, adaptive.frac, trace, tol, pos, thr = cor.thr, early_flag = T)
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



    #md.run <- md.shrink(Y, B0, k, max.iter, L1, L2, L4, right.shur, adaptive.iter, adaptive.frac, trace, tol)
    #while (md.run$flag == "Shrink"){
    #  L4 <- L4/L4_adaptive
    #  print(L4)
    #  md.run <- md.shrink(Y, B0, k, max.iter, L1, L2, L4, right.shur, adaptive.iter, adaptive.frac, trace, tol)
    #}#while

  #wrap around the output
  #######################################################
  B <- md.run$B
  Z <- md.run$Z
  Zraw <- md.run$Zraw
  Zproject=Z%*%solve(crossprod(Z)+L2*diag(k))

  if (to_drop){

    if (verbose){
      message("drop")
    }#if verbose

    #eliminate the one with smaller variance
    cor.res <- cor(t(B))
    LV.var <- VarianceExplained(Y, Z, B, option = "simple", normalize = F)
    drop.index <- which(cor.res == max(cor.res[upper.tri(cor.res)]) & upper.tri(cor.res), arr.ind = T)
    if (nrow(drop.index)!=1){
      warning(nrow(drop.index))
    }#if
    LV.to_drop <- drop.index[which.min(LV.var[drop.index])]

    #keep a record
    B0 <- B
    Z0 <- Z

    #exclue LV.to_drop
    B <- B[-LV.to_drop,]
    Z <- Z[,-LV.to_drop]
    Zraw <- Zraw[,-LV.to_drop]
    Zproject <- Zproject[,-LV.to_drop]
    k <- k-1

  }#to_drop

  #assign names
  rownames(B)=colnames(Z)=paste("LV",1:k)
  rownames(Z) <- rownames(Y)
  colnames(B) <- colnames(Y)

  if (save.complete){
    if (to_drop){
      return(list(B=B, Z=Z, Zraw=Zraw, Zproject=Zproject,L1=L1, L2=L2, L4 = L4_pointer, k = k, shur0 = shur0, right.shur = right.shur, Y= Y, L= L, LV.to_drop = LV.to_drop, B0 = B0, Z0 = Z0))
    }else{
      return(list(B=B, Z=Z, Zraw=Zraw, Zproject=Zproject,L1=L1, L2=L2, L4 = L4_pointer, k = k, shur0 = shur0, right.shur = right.shur, Y= Y, L= L))
    }#else

  }else{
    return(list(B=B, Z=Z, Zraw=Zraw, Zproject=Zproject,L1=L1, L2=L2, L4 = L4_pointer, k = k, shur0 = shur0, right.shur = right.shur))
  }#else


}#manifoldDecomp_adaptive











#tSNE or Umap plot
#' Generate 2D projection plot
#' @import ggpubr
#' @details This function requires `uwot` and `Rtsne`. Make sure they are installed.
#' @export
Cluster_2Dplot <- function(LVs, cluster.label = NULL, cluster.option = "kmeans", cluster.label.scale = F, num.of.cluster = NULL,option = "umap", random.seed = 1, pt.size = 2, verbose = T, umap.opt = 30, tsne.opt = 30){

  message(paste0("Generate ", option, " plot."))
  set.seed(random.seed)

  kmeans.input <- t(LVs)
  if (is.null(cluster.label)){

    if (cluster.label.scale){
      kmeans.input <- scale(t(LVs))
    }#if

    if (is.null(num.of.cluster)){
      stop("num.of.cluster is unknown.")
    }#if

    if (cluster.option == "kmeans"){
      km.res <- kmeans(kmeans.input, centers = num.of.cluster, nstart = 25)
      cluster.label <- as.character(km.res$cluster)
    }else if (cluster.option == "leiden"){

    }#else if

  }#if is.null(cluster.label)


  #visualization
  if (option == "umap"){

    if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("uwot is required for this function but is not installed. Please install it.")
  }#if

    #umap.proj <- uwot::umap(kmeans.input, n_neighbors = 30, learning_rate = 0.5, init = "random", verbose = F, n_threads = 4, pca = NULL, n_epochs = 1000)
    umap.proj <- uwot::umap(kmeans.input, n_neighbors = umap.opt, learning_rate = 0.5, init = "random", verbose = F, n_threads = 4, pca = NULL, n_epochs = 1000)
    dat <- data.frame(x = umap.proj[,1], y = umap.proj[,2], Cluster = as.character(cluster.label))

  }else if (option == "tsne"){

    if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Rtsne is required for this function but is not installed. Please install it.")
  }#if

    #rtsne.proj <- Rtsne::Rtsne(kmeans.input, perplexity = 30, verbose = F, pca = F)
    rtsne.proj <- Rtsne::Rtsne(kmeans.input, perplexity = tsne.opt, verbose = F, pca = F)
    dat <- data.frame(x = rtsne.proj$Y[,1], y = rtsne.proj$Y[,2], Cluster = as.character(cluster.label))
  }else if (option == "pca"){
    svd.res  <- svd(t(scale(kmeans.input)), nu = 2, nv = 2)
    dat <- data.frame(x = svd.res$v[,1], y = svd.res$v[,2], Cluster = as.character(cluster.label))
  }#else if

  rownames(dat) <- rownames(kmeans.input)
  if (verbose){
    ggplot(dat,aes(x=x,y=y, color = Cluster))+geom_point(size = pt.size)+theme_pubr(base_size = 20)
  }else{
    return(dat)
  }#else

}#Cluster_2Dplot






