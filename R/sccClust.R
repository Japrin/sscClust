#' Single cell RNA-Seq data extracted from a publication by Pollen et al.
#'
#' @name sce.Pollen
#' @docType data
#' @source \url{https://www.nature.com/articles/nbt.2967}
#'
#' Expression data contained in a \code{SingleCellExperiment} object
NULL

#' dispaly message with time stamp
#' @param msg characters; message to display
loginfo <- function(msg) {
  timestamp <- sprintf("%s", Sys.time())
  msg <- paste0("[",timestamp, "] ", msg,"\n")
  cat(msg)
}

#' Generate color set automatically
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @param n integer; number of colors needed
#' @param name character; Palette name
#' @return an vector contains the color codes
auto.colSet <- function(n=2,name="Set1"){
  requireNamespace("RColorBrewer",quietly = T)
  if(n<=8){
    ret <- RColorBrewer::brewer.pal(max(n,3),name)[seq_len(n)]
  }else{
    ret <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(n)
  }
  return(ret)
}

#' Determine point size automatically
#' @param n number of points to plot
#' @return points' cex
auto.point.size <- function(n){
  if(n<=100){
    return(1.2)
  }else if(n>=5000){
    return(0.6)
  }else{
    return(-0.6*n/4900+1.212002)
  }
}

#' Plot gene expression on tSNE map
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot ggsave scale_colour_gradientn geom_point facet_wrap theme_bw
#' @importFrom data.table melt
#' @importFrom utils head
#' @param Y matrix or data.frame; Gene expression data, rownames shoud be gene id, colnames
#' should be sample id
#' @param dat.map data.frame; tSNE map, must be two columns data.frame and rownames should be sample id
#' @param gene.to.show character; gene id to be showd on the tSNE map
#' @param out.prefix character; output prefix (default: NULL)
#' @param p.ncol integer; number of columns in the plot's layout (default: 3)
#' @param width numeric; width of the plot (default: 9)
#' @param height numeric; height of the plot (default: 8)
#' @details For genes contained in both `Y` and `gene.to.show`, show their expression on the tSNE
#' map provided as `dat.map`. One point in the map represent a cell; cells with higher expression
#' also have darker color.
#' @return a ggplot object
ggGeneOnTSNE <- function(Y,dat.map,gene.to.show,out.prefix=NULL,p.ncol=3,width=9,height=8){
  #suppressPackageStartupMessages(require("data.table"))
  #requireNamespace("ggplot2",quietly = T)
  #requireNamespace("RColorBrewer",quietly = T)
  if(!is.null(out.prefix)){
    dir.create(sprintf("%s.perGene.tSNE",out.prefix),showWarnings = F,recursive = T)
  }
  f.g <- gene.to.show %in% rownames(Y)
  if(sum(!f.g)>0){
    warning(sprintf("Some genes not in the expression data: \n"))
    print(head(gene.to.show[!f.g]))
  }
  gene.to.show <- gene.to.show[f.g]

  dat.plot <- data.frame(sample=rownames(dat.map),stringsAsFactors = F)
  dat.plot <- cbind(dat.plot,dat.map,t(Y[gene.to.show,dat.plot$sample,drop=F]))
  colnames(dat.plot) <- c("sample","Dim1","Dim2",names(gene.to.show))
  dat.plot.melt <- data.table::melt(dat.plot,id.vars = c("sample","Dim1","Dim2"))
  dat.plot.melt <- dat.plot.melt[order(dat.plot.melt$value,decreasing = F),]
  npts <- nrow(dat.plot.melt)
  p <- ggplot2::ggplot(dat.plot.melt,aes(Dim1,Dim2)) +
          geom_point(aes(colour=value),size=auto.point.size(npts)*1.1) +
          scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd")) +
          facet_wrap(~variable, ncol = p.ncol) +
          theme_bw()
  if(!is.null(out.prefix)){
    ggplot2::ggsave(sprintf("%s.geneOntSNE.pdf",out.prefix),width = width,height = height)
  }
  return(p)
}

#' Find the knee point of the scree plot
#'
#' @param pcs principal component values sorted decreasingly
#' @details Given sorted decreasingly PCs, find the knee point which have the largest distance to
#' the line defined by the first point and the last point in the scree plot
#' @return index of the knee plot
findKneePoint <- function(pcs)
{
  npts <- length(pcs)
  if(npts<=3){
    return(npts)
  }else{
    P1 <- c(1,pcs[1])
    P2 <- c(npts,pcs[npts])
    v1 <- P1 - P2
    dd <- sapply(2:(npts-1),function(i){
      Pi <- c(i, pcs[i])
      v2 <- Pi - P1
      m <- cbind(v1,v2)
      d <- abs(det(m))/sqrt(sum(v1*v1))
    })
    return(which.max(dd))
  }
}

####### classification functions

#' Wraper for running random forest classifier
#'
#' @importFrom varSelRF varSelRF
#' @importFrom stats predict
#' @param xdata data frame or matrix; data used for training, with sample id in rows and variables in columns
#' @param xlabel factor; classification label of the samples, with length equal to the number of rows in xdata
#' @param ydata data frame or matrix; data to be predicted the label, same format as xdata
#' @param do.norm logical; whether perform Z score normalization on data
#' @return List with the following elements:
#' \item{ylabel}{ppredicted labels of the samples in ydata}
#' \item{rfsel}{trained model; output of varSelRF()}
run.RF <- function(xdata, xlabel, ydata, do.norm=F)
{
  #require("varSelRF")
  #require("randomForest")
  f.g <- intersect(colnames(xdata),colnames(ydata))
  xdata <- xdata[,f.g,drop=F]
  ydata <- ydata[,f.g,drop=F]
  ### normalization
  if(do.norm){
    xdata <- scale(xdata,center = T,scale = T)
    ydata <- scale(ydata,center = T,scale = T)
  }
  ### random forest
  rfsel <- varSelRF::varSelRF(xdata, xlabel,
                    ntree = 500, ntreeIterat = 200,
                    whole.range = FALSE,keep.forest = T)
  #rfsel$selected.vars %>% str %>% print
  #rfsel$initialImportances %>% head %>% print
  #rfsel$rf.model$confusion %>% print
  ylabel <- predict(rfsel$rf.model, newdata = ydata[,rfsel$selected.vars])
  names(ylabel) <- rownames(ydata)
  return(list("ylabel"=ylabel,"rfsel"=rfsel))
}

#' Wraper for running random forest classifier
#' @importFrom class knn
#' @param xdata data frame or matrix; data used for training, with sample id in rows and variables in columns
#' @param xlabel factor; classification label of the samples, with length equal to the number of rows in xdata
#' @param ydata data frame or matrix; data to be predicted the label, same format as xdata
#' @param k parameter k of function knn() (default: 1)
#' @return List with the following elements:
#' \item{ylabel}{ppredicted labels of the samples in ydata}
run.KNN <- function(xdata,xlabel,ydata,k=1)
{
  #require("class")
  f.g <- intersect(colnames(xdata),colnames(ydata))
  ylabel <- class::knn(xdata[,f.g,drop=F], ydata[,f.g,drop=F], as.factor(xlabel), k = k, l = 0, prob = FALSE, use.all = TRUE)
  names(ylabel) <- rownames(ydata)
  return(list("ylabel"=ylabel))
}

#' Wraper for running Rtsne
#' @importFrom Rtsne Rtsne
#' @param idata matrix; expression data with sample id in rows and variables in columns
#' @param tSNE.usePCA whether perform PCA before tSNE (default: T)
#' @param tSNE.perplexity perplexity parameter of tSNE (default: 30)
#' @return If successful same as the return value of Rtsne(); otherwise NULL
run.tSNE <- function(idata,tSNE.usePCA=T,tSNE.perplexity=30){
  ret <- NULL
  tryCatch({
    ret <- Rtsne::Rtsne(idata, pca = tSNE.usePCA, perplexity = tSNE.perplexity)$Y
  },error=function(e){
      #cat("Perplexity is too large; try to use smaller perplexity 5\n")
    })
  if(is.null(ret)){
    tryCatch({
      ret <- Rtsne::Rtsne(idata, pca = tSNE.usePCA, perplexity = 5)$Y
    },error=function(e){ print("Error occur when using perplexity 5"); e })
  }
  return(ret)
}

#' Wraper for silhouette()
#' @importFrom cluster silhouette
#' @importFrom graphics plot
#' @importFrom stats dist
#' @param obj object of \code{SingleCellExperiment}
#' @param cluster.label character; which column of colData of obj to used as cluster label.
#' @param reducedDim.name character; which reducedDim to use. (default: "iCor.tsne")
#' @param do.plot logical; whether plot
#' @param ... Arguments to be passed to plot()
#' @return an object, sil, of class \code{silhouette}
#' @export
ssc.plot.silhouette <- function(obj,cluster.label,reducedDim.name="iCor.tsne",do.plot=T, ...){
  requireNamespace("cluster")
  dist.obj <- dist(reducedDim(obj,reducedDim.name))
  sil <- cluster::silhouette(as.numeric(as.factor(colData(obj)[,cluster.label])),dist.obj)
  if(do.plot){
    plot(sil, ...)
  }
  return(sil)
}


#### ===================================

#' Build an SingleCellExperiment object
#'
#' Build an SingleCellExperiment object from a matrix or data frame
#'
#' @importFrom SingleCellExperiment SingleCellExperiment rowData
#' @param x matrix/data.frame or SingleCellExperiment; input expression data
#' @param display.name a vector, should be human readable gene name
#' @details if x is an object of SingleCellExperiment, just clear the metadata;
#' if x is matrix/data.frame, convert it to an object of SingleCellExperiment.
#' Also a vector `display.name` can be provided, which would be used in some plots,
#' such as geneOnTSNE, heatmap. The row names of SingleCellExperiment object usually be gene id
#' (e.g. entrez ID, Ensemble ID), the `display.name` should be human readable gene name (
#' e.g. HGNC gene symbol). If `display.name` is NULL (default), the row names of SingleCellExperiment object
#' would be used.
#' @return an object of \code{SingleCellExperiment} class
#' @export
ssc.build <- function(x,display.name=NULL)
{
  obj <- NULL
  if(class(x)=="SingleCellExperiment")
  {
    obj <- x
    metadata(obj)$ssc <- list()
  }else if(class(x) %in% c("matrix","data.frame")){
    obj <- SingleCellExperiment(assays = list(exprs = as.matrix(x)))
    if(!is.null(display.name)){
      f.na <- is.na(display.name)
      display.name[f.na] <- row.names(obj)[f.na]
      rowData(obj)[,"display.name"] <- display.name
    }else{
      rowData(obj)[,"display.name"] <- row.names(obj)
    }
  }
  return(obj)
}

#' Identify variable genes
#'
#' Identify variable genes which will be used in downstream analysis. Multiple methods are available.
#'
#' @importFrom stats sd
#' @param obj object of SingleCellExperiment
#' @param method method to be used, can be one of "sd" (default),
#' @param sd.n top number of genes (default 1500)
#' @param assay.name which assay to be used (default "exprs")
#' @details Method "sd", calculate the standard deviation of each gene and sort decreasingly, the top `sd.n` genes are the
#' variable genes.
#' @return an object of \code{SingleCellExperiment} class
#' @export
ssc.variableGene <- function(obj,method="sd",sd.n=1500,assay.name="exprs")
{
  if(method=="sd")
  {
    row.sd <- apply(assay(obj,assay.name),1,sd)
    metadata(obj)$ssc[["variable.gene"]][["sd"]] <- names(head(row.sd[order(row.sd,decreasing = T)],n=sd.n))
  }else if(method=="mean.sd")
  {
    row.sd <- apply(assay(obj,assay.name),1,sd)
    row.mean <- apply(assay(obj,assay.name),1,mean)
    info.gene.sd <- row.sd[row.sd>1 & row.mean>1]
    info.gene.sd <- sort(info.gene.sd,decreasing = T)
    metadata(obj)$ssc[["variable.gene"]][["mean.sd"]] <- names(head(info.gene.sd,n=sd.n))
  }
  return(obj)
}

#' Convert gene id to display name (gene symbol)
#' @importFrom SingleCellExperiment rowData
#' @param obj object of \code{SingleCellExperiment} class
#' @param ids character; gene ids
#' @return If successfull vector contains the disaply name; otherwise NULL
#' @export
ssc.id2displayName <- function(obj,ids)
{
  ret <- NULL
  if("display.name" %in% colnames(rowData(obj)))
  {
    lookup.table <- structure(rowData(obj)[,"display.name"],
                              names=rownames(obj))
    ret <- lookup.table[ids]
  }
  return(ret)
}

#' Convert display name (gene symbol) to gene id
#' @importFrom SingleCellExperiment rowData
#' @param obj object of \code{SingleCellExperiment} class
#' @param display.name character; disaply name
#' @return If successfull vector contains the gene id; otherwise NULL
#' @export
ssc.displayName2id <- function(obj,display.name)
{
  ret <- NULL
  if("display.name" %in% colnames(rowData(obj)))
  {
    lookup.table <- structure(rowData(obj)[,"display.name"],
                              names=rownames(obj))
    #ret <- lookup.table[ids]
    ret <- lookup.table[lookup.table %in% display.name]
    ret <- structure(names(ret),names=lookup.table[names(ret)])
  }
  return(ret)
}

#' Reduce dimension by various methods
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom stats cor prcomp
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param method character; method to be used for dimension reduction, should be one of (pca, tsne, iCor). (default: "iCor")
#' @param method.vgene character; method to identify variable genes. (default: sd)
#' @param pca.npc integer; number of pc be used. Only for reduction method "pca". (default: NULL)
#' @param tSNE.usePCA logical; whether use PCA before tSNE. Only for reduction method "tsne". (default: T)
#' @param tSNE.perplexity logical; perplexity parameter. Used in all Rtsne() calling. (default: 30)
#' @param autoTSNE logical; Wheter generate automatically a tSNE map when reduction method is "pca" or "iCor". (default: T)
#' @param iCor.niter integer; number of iteration of calculating the correlation. Used in reduction method "iCor". (default: 1)
#' @param iCor.method character; method to calculate correlation between samples,
#' should be one of "spearman" and "pearson". (default "spearman")
#' @details If the reduction method is "pca", the function will call prcomp() and estimate the number of top PC should be used
#' in downstream analysis using and "elbow" based method, then the samples coordinates in the space spaned by the top PC would
#' stored in the reducedDim slot of the return value with the reducedDimName "pca".If autoTSNE is `true`, a tSNE map based on
#' the top PC will generated and stored in the reducedDim slot with reducedDimName "pca.tsne".
#'
#' If the reduction method is "tsne", a tSNE map will generated and stored in the reducedDim slot with reducedDimName "tsne".
#' If tSNE.usePCA is `true`, maximum 50 top PC will be used, but no estimation based on "elbow" method performed.
#'
#' If the reduction method is "iCor", the function will calculate correlation between samples for iCor.niter times.
#' If autoTSNE is `true`, it also generate automatically a tSNE map based on the correlations. These data would be stored in the
#' reducedDim slot with reducedDimName "iCor" and "iCor.tsne" respectively.
#' @return an object of \code{SingleCellExperiment} class with reduced data added to the reducedDim slot.
#' @export
ssc.reduceDim <- function(obj,assay.name="exprs",
                          method="iCor",
                          method.vgene="sd",
                          pca.npc=NULL,
                          tSNE.usePCA=T,
                          tSNE.perplexity=30,
                          autoTSNE=T,
                          iCor.niter=1,iCor.method="spearman")
{
  row.sd <- apply(assay(obj,assay.name),1,sd)
  col.sd <- apply(assay(obj,assay.name),2,sd)
  col.zero <- which(col.sd==0)
  if(length(col.zero>0)){
    warning(sprintf("expression data contains colum(s) with sd equal to 0:\n%s\n",
                    paste(head(colnames(obj)[col.zero]),collapse = ",")))
  }
  obj <- obj[row.sd>0,]
  if(!method.vgene %in% names(metadata(obj)$ssc[["variable.gene"]])){
    stop(sprintf("No variable genes identified by method %s !",method.vgene))
  }
  vgene <- metadata(obj)$ssc[["variable.gene"]][[method.vgene]]
  if(method=="pca"){
    pca.res <- prcomp(t(assay(obj[vgene,],assay.name)))
    ### find elbow point and get number of components to be used
    pca.res$eigenv.prop <- pca.res$sdev/sum(pca.res$sdev)
    pca.res$eigengap <- sapply(seq_len(length(pca.res$eigenv.prop)-1),function(i){ pca.res$eigenv.prop[i]-pca.res$eigenv.prop[i+1] })
    ### method 1
    ######pca.res$kneePts <- which(pca.res$eigengap<1e-4)[1]
    ### method 2 (max distance to the line defined by the first and last point in the scree plot)
    pca.res$kneePts <- findKneePoint(head(pca.res$eigenv.prop,n=50))
    if(is.null(pca.npc) && !is.na(pca.res$kneePts)){ pca.npc <- pca.res$kneePts
    }else{ pca.npc <- 30 }
    pca.npc <- min(pca.npc,ncol(pca.res$x))
    pca.res$npc <- pca.npc
    ### save to object
    metadata(obj)$ssc$pca.res <- pca.res
    proj_data <- pca.res$x[,1:pca.npc,drop=F]
    if(autoTSNE){ reducedDim(obj,"pca.tsne") <- run.tSNE(proj_data,tSNE.usePCA=F,tSNE.perplexity) }
  }else if(method=="tsne"){
    proj_data <- run.tSNE(t(assay(obj[vgene,],assay.name)),tSNE.usePCA,tSNE.perplexity)
  }else if(method=="iCor"){
    proj_data <- assay(obj[vgene,],assay.name)
    while(iCor.niter>0){
      proj_data <- cor(proj_data,method=iCor.method)
      iCor.niter <- iCor.niter-1
    }
    if(autoTSNE) { reducedDim(obj,"iCor.tsne") <- run.tSNE(proj_data,tSNE.usePCA=F,tSNE.perplexity) }
  }
  reducedDim(obj,method) <- proj_data
  return(obj)
}


#' Perform clustering using reduced data
#' @importFrom factoextra eclust
#' @importFrom stats kmeans
#' @importFrom ADPclust adpclust
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param method.reduction character; which dimention reduction method to be used, should be one of
#' "iCor", "pca" and "none". (default: "iCor")
#' @param method character; clustering method to be used, should be one of "kmeans", "hclust" an "adpclust". (default: "kmeans")
#' @param k.batch integer; number of clusters to be evaluated. (default: 2:6)
#' @param method.vgene character; variable gene identification method used. (default: "sd")
#' @details If no dimension reduction performed or method is "none", expression data of variable genes,
#' which can be speficed by method.vgene, will be used for clustering. Otherwise, the reduced data specified by
#' method.reduction will be used. The cluster label will stored in the colData of the object
#' of \code{singleCellExperiment} class, with colname in the format of \{method.reduction\}.\{method\}k\{k\}
#' where \{k\} get value(s) from k.batch.
#' @return an object of \code{SingleCellExperiment} class with cluster labels added.
#' @export
ssc.clust <- function(obj, assay.name="exprs", method.reduction="iCor",
                      method="kmeans", k.batch=2:6,
                      method.vgene="sd")
{
  clust.res <- NULL
  res.list <- list()
  ### check transformed data
  if(method.reduction=="none" || (!method.reduction %in% reducedDimNames(obj)) ){
    warning(sprintf("The dimention reduction should be performed!"))
    vgene <- metadata(obj)$ssc[["variable.gene"]][[method.vgene]]
    dat.transformed <- t(assay(obj[vgene,],assay.name))
  }else{
    dat.transformed <- reducedDim(obj,method.reduction)
  }
  ### check method

  for(k in k.batch){
    if(method=="kmeans"){
      clust.res <- kmeans(dat.transformed,k,iter.max=1000,nstart=50)
      colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%d",clust.res$cluster)
    }else if(method=="hclust"){
      clust.res <- factoextra::eclust(dat.transformed, "hclust", k = k, method = "complete", graph = FALSE)
      colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%d",clust.res$cluster)
    }else if(method=="adpclust"){
      clust.res <- ADPclust::adpclust(dat.transformed,nclust = k.batch)
      k <- "auto"
      colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%d",clust.res$clusters)
    }
    res.list[[as.character(k)]] <- clust.res
    if(method=="adpclust"){ break }
  }
  metadata(obj)$ssc$clust.res[[method]] <- res.list
  return(obj)
}

#' Clustering with subsampling and classification
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom stats cor
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param frac numeric; subsample to frac of original samples. (default: 0.4)
#' @param method.vgene character; variable gene identification method used. (default: "sd")
#' @param method.reduction character; which dimention reduction method to be used, should be one of
#' "iCor", "pca" and "none". (default: "iCor")
#' @param method.clust character; clustering method to be used, should be one of "kmeans" and "hclust". (default: "kmeans")
#' @param method.classify character; method used for classification, one of "knn" and "RF". (default: "knn")
#' @param pca.npc integer; number of pc be used. Only for reduction method "pca". (default: NULL)
#' @param iCor.niter integer; number of iteration of calculating the correlation. Used in reduction method "iCor". (default: 1)
#' @param use.proj logical; whether use the projected data for classification. (default: T)
#' @param k.batch integer; number of clusters to be evaluated. (default: 2:6)
#' @details The function first subsmaple the samples to the specified fraction (such 40%), and perform clustering. The clustering will
#' make labels for the subsampled samples. Using the labels, original data or projected data via the method specified in
#' "method.reduction" will be used for trainning a classifier. Then the classifier will predict the labels of the samples not subsampled,
#' using original data or projected data dependent on the option use.proj. The final cluster labels combining
#' that of bath sampled and unsampled samples, will stored in the colData of the object of \code{singleCellExperiment} class,
#' with colname in the format of \{method.reduction\}.\{method\}k\{k\} where \{k\} get value(s) from k.batch.
#' @return an object of \code{SingleCellExperiment} class with cluster labels added.
#' @export
ssc.clustSubsamplingClassification <- function(obj, assay.name="exprs",
                                               frac=0.4,
                                               method.vgene="sd",
                                               method.reduction="iCor",
                                               method.clust="kmeans",
                                               method.classify="knn",
                                               pca.npc=NULL,
                                               iCor.niter=1,
                                               use.proj=T,
                                               k.batch=2:6)
{
  #### subsampling
  n.sub <- floor(ncol(obj)*frac)
  if(n.sub<10){
    stop(sprintf("too few samples, n.sub: %d\n",n.sub))
  }
  cids <- sample(ncol(obj),n.sub)
  obj.train <- obj[,cids]
  obj.pred <- obj[,-cids]
  #### because the samples changed, obj.train/obj.pred may contain genes with zero sd
  f.sd0 <- (apply(assay(obj.train,assay.name),1,sd)==0 | apply(assay(obj.pred,assay.name),1,sd)==0)
  obj.train <- obj.train[!f.sd0,]
  obj.pred <- obj.pred[!f.sd0,]

  reducedDims(obj.train) <- SimpleList()
  #### reduction
  if(method.reduction=="none"){
    warning(sprintf("Reducing the dimensions first is recommended"))
  }else if(method.reduction %in% c("pca","iCor")){
    obj.train <- ssc.reduceDim(obj.train,method=method.reduction,
                              method.vgene=method.vgene,
                              pca.npc=pca.npc,autoTSNE = T,
                              iCor.niter=iCor.niter,iCor.method="spearman")
    if(!all(rownames(obj.train)==rownames(obj.pred))){
      stop("The genes of obj.train and obj.pred are different!")
    }
    vgene <- metadata(obj.pred)$ssc$variable.gene[[method.vgene]]

    if(method.reduction=="iCor"){
      #### for (tsne) visualization
      A <- cor((assay(obj.pred[vgene,],assay.name)),assay(obj.train[vgene,],assay.name),method = "spearman")
      S <- cor(assay(obj.train[vgene,],assay.name),method = "spearman")
      #S <- reducedDim(obj.train,sprintf("%s",method.reduction))
      dat.map.train <- reducedDim(obj.train,sprintf("%s.tsne",method.reduction))
      dat.map.pred <- A %*% solve(S) %*% dat.map.train
      reducedDim(obj.pred,sprintf("%s.tsne",method.reduction)) <- dat.map.pred
      dat.map.obj <- rbind(dat.map.train,dat.map.pred)
      reducedDim(obj,sprintf("%s.tsne",method.reduction)) <- dat.map.obj[colnames(obj),,drop=F]
      #### for prediction
      metadata(obj.pred)$ssc$data.transformed <- A
    }else if(method.reduction=="pca"){
      #### for (tsne) visualization
      #### for prediction
      metadata(obj.pred)$ssc$data.transformed <- t(assay(obj.pred[vgene,],assay.name)) %*% metadata(obj.train)$ssc$pca.res$rotation
      metadata(obj.pred)$ssc$data.transformed <- metadata(obj.pred)$ssc$data.transformed[,1:metadata(obj.train)$ssc$pca.res$npc,drop=F]
    }
  }else{
    stop(sprintf("unsupported dimension reduction method: %s\n",method.reduction))
  }
  #### clustering
  obj.train <- ssc.clust(obj.train, method.reduction=method.reduction,
                         method=method.clust, k.batch=k.batch)
  colData(obj)[,"isTrainSet"] <- F
  colData(obj)[colnames(obj.train),"isTrainSet"] <- T
  #### data for classification
  if(use.proj && method.reduction %in% c("pca","iCor")){
    data.train <- reducedDim(obj.train,sprintf("%s",method.reduction))
    data.pred <- metadata(obj.pred)$ssc$data.transformed
  }else{
    data.train <- t(assay(obj.train,assay.name))
    data.pred <- t(assay(obj.pred,assay.name))
  }
  #### classificaton
  ## option 1: determin k.best
  #cls.label <- colData(obj.train)[,sprintf("%s.%s.k%s",method.reduction,method.clust,k.best)]
  ## option 2: loop all k.batch
  ## parallel
  for(kk in k.batch){
    cls.label <- colData(obj.train)[,sprintf("%s.%s.k%s",method.reduction,method.clust,kk)]
    names(cls.label) <- colnames(obj.train)
    if(method.classify=="RF")
    {
      res.pred <- run.RF(data.train, as.factor(cls.label), data.pred,do.norm = T)
    }else if(method.classify=="knn")
    {
      res.pred <- run.KNN(data.train,as.factor(cls.label), data.pred,k=1)
    }else if(method.classify=="svm")
    {
    }
    pred.label <- structure(as.character(res.pred$ylabel), names=names(res.pred$ylabel))
    ret.label <- c(cls.label,pred.label)
    colData(obj)[names(ret.label),sprintf("%s.%s.k%s",method.reduction,method.clust,kk)] <- ret.label
  }
  return(obj)
}

#' Wrapper for running all the pipeline
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param method.vgene character; variable gene identification method used. (default: "sd")
#' @param sd.n top number of genes (default 1500)
#' @param method.reduction character; which dimention reduction method to be used, should be one of
#' "iCor", "pca" and "none". (default: "iCor")
#' @param method.clust character; clustering method to be used, should be one of "kmeans" and "hclust". (default: "kmeans")
#' @param method.classify character; method used for classification, one of "knn" and "RF". (default: "knn")
#' @param pca.npc integer; number of pc be used. Only for reduction method "pca". (default: NULL)
#' @param iCor.niter integer; number of iteration of calculating the correlation. Used in reduction method "iCor". (default: 1)
#' @param subsampling logical; whether cluster using the subsampling->cluster->classification method. (default: F)
#' @param sub.frac numeric; subsample to frac of original samples. (default: 0.4)
#' @param sub.use.proj logical; whether use the projected data for classification. (default: T)
#' @param k.batch integer; number of clusters to be evaluated. (default: 2:6)
#' @details run the pipeline of variable gene identification, dimension reduction, clustering.
#' @seealso \code{\link{ssc.variableGene}} for variable genes' identification, \code{\link{ssc.reduceDim}}
#' for dimension reduction, \code{\link{ssc.clust}} for clustering using all data
#' and \code{\link{ssc.clustSubsamplingClassification}} for clustering with subsampling.
#' @return an object of \code{SingleCellExperiment} class with cluster labels added.
#' @export
ssc.run <- function(obj, assay.name="exprs",
                    method.vgene="sd",
                    sd.n=1500,
                    method.reduction="iCor",
                    method.clust="kmeans",
                    method.classify="knn",
                    pca.npc=NULL,
                    iCor.niter=1,
                    subsampling=F,
                    sub.frac=0.4,
                    sub.use.proj=T,
                    k.batch=2:6)
{
  ### some checking
  if(is.null(colnames(obj))){
    stop("colnames of obj is NULL!!!")
  }
  if(is.null(rownames(obj))){
    stop("rownames of obj is NULL!!!")
  }

  obj <- ssc.variableGene(obj,method=method.vgene,sd.n=sd.n,assay.name=assay.name)
  if(!subsampling){
    obj <- ssc.reduceDim(obj,assay.name=assay.name,
                              method=method.reduction,
                              pca.npc = pca.npc,
                              iCor.niter = iCor.niter,
                              method.vgene=method.vgene)
    obj <- ssc.clust(obj, assay.name=assay.name, method.reduction=method.reduction,
                          method=method.clust, k.batch=k.batch,
                          method.vgene=method.vgene)
  }else{
    obj <- ssc.clustSubsamplingClassification(obj, assay.name=assay.name,
                                                   frac=sub.frac,
                                                   method.vgene=method.vgene,
                                                   method.reduction=method.reduction,
                                                   method.clust=method.clust,
                                                   method.classify=method.classify,
                                                   pca.npc=pca.npc,
                                                   iCor.niter=iCor.niter,
                                                   use.proj=sub.use.proj,
                                                   k.batch=k.batch)
  }
  return(obj)
}

#' Plot on tSNE map
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param gene, character; genes to be showed. (default: NULL)
#' @param columns character; columns in colData(obj) to be showd. (default: NULL)
#' @param colSet list; mapping iterms in the names to colors in the values. (default: list())
#' @param reduced.name character; names in the reducedDimNames. (default: "iCor.tsne")
#' @param reduced.dim integer; which dimensions of the reduced data to be used. (default: c(1,2))
#' @param out.prefix character; output prefix. (default: NULL)
#' @param p.ncol integer; number of columns in the figure layout. (default: 3)
#' @param width numeric; width of the plot, used for geneOnTSNE. (default: NA)
#' @param height numeric; height of the plot, used for geneOnTSNE. (default: NA)
#' @param base_aspect_ratio numeric; base_aspect_ratio, used for plotting metadata. (default 1.1)
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual theme_bw aes_string guides guide_legend
#' @importFrom cowplot save_plot plot_grid
#' @importFrom utils read.table
#' @details If `gene` is not NULL, expression of the specified genes will be plot on the tSNE map; if columns in not
#' NULL, colData of obj with names in `columns` will be plot on the tSNE map. The tSNE map used is specified by option
#' `reduced.name` and `reduced.dim`. Both `gene` and `columns` can be non-NULL. For list `colSet`, each element define
#' a color mapping for the responding iterm in the `column`; if not specifed, automatically generated color mapping will
#' be used.
#' @export
ssc.plot.tsne <- function(obj, assay.name="exprs", gene=NULL, columns=NULL, colSet=list(),
                          reduced.name="iCor.tsne",reduced.dim=c(1,2),
                          out.prefix=NULL,p.ncol=3,width=NA,height=NA,base_aspect_ratio=1.1)
{
  #requireNamespace("ggplot2")
  #requireNamespace("cowplot")
  if(length(reduced.dim)!=2){ stop(sprintf("Wrong parameter, reduced.dim!!"))}
  if(!is.null(columns))
  {
    if(all(columns %in% colnames(colData(obj))))
    {
      if(is.list(colSet)){
        dat.map <- reducedDim(obj,reduced.name)
        multi.p <- lapply(columns,function(cc){
          if(is.null(colSet[[cc]])){
            cc.values <- sort(unique(colData(obj)[,cc]))
            colSet[[cc]] <- structure(auto.colSet(length(cc.values),name = "Paired"),
                                      names=cc.values)
          }
          dat.plot <- data.frame(sample=rownames(dat.map),stringsAsFactors = F)
          dat.plot <- as.data.frame(cbind(dat.plot,dat.map[,reduced.dim],colData(obj)[,cc,drop=F]))
          colnames(dat.plot) <- c("sample","Dim1","Dim2",cc)
          npts <- nrow(dat.plot)
          p <- ggplot2::ggplot(dat.plot,aes(Dim1,Dim2)) +
            geom_point(aes_string(colour=cc),size=auto.point.size(npts)*1.1) +
            scale_colour_manual(values = colSet[[cc]]) +
            theme_bw() +
            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=4)))
          return(p)
        })
        pp <- cowplot::plot_grid(plotlist=multi.p,ncol = if(length(columns)>1) 2 else 1,align = "hv")
        if(!is.null(out.prefix)){
          cowplot::save_plot(sprintf("%s.columnsOntSNE.pdf",out.prefix),pp,
                             ncol = if(length(columns)>1) 2 else 1,
                             base_aspect_ratio=base_aspect_ratio)
        }else{
          print(pp)
        }
      }else{
        stop(sprintf("invalidate parameter: colSet. Please check that!"))
      }
    }else{
      warning(sprintf("some columns not in the data. Not plot be produced!"))
    }
  }
  if(!is.null(gene)){
    if(length(gene)==1 && file.exists(gene)){ gene <- read.table(gene,header = T)[,1] }
    if(all(!(gene %in% rownames(obj)))){
      ### gene symbol?
      gene <- ssc.displayName2id(obj,display.name = gene)
    }
    p <- ggGeneOnTSNE(assay(obj,assay.name),
                     reducedDim(obj,reduced.name)[,reduced.dim],
                     gene,out.prefix,p.ncol=p.ncol,width=width,height=height)
    if(is.null(out.prefix)){ print(p) }
  }
}

#' Plot pca result, such as scree plot.
#' @param obj object of \code{singleCellExperiment} class
#' @param out.prefix character; output prefix. (default: NULL)
#' @param p.ncol integer; number of columns in the figure layout. (default: 2)
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual theme_bw
#' @export
ssc.plot.pca <- function(obj, out.prefix=NULL,p.ncol=2){
  requireNamespace("ggplot2")
  eigenv <- metadata(obj)$ssc$pca.res$eigenv.prop
  dat.plot.eigenv <- data.frame(PC=seq_along(eigenv),
                                eigenv=eigenv,
                                isKneePts=as.character(seq_along(eigenv)==metadata(obj)$ssc$pca.res$kneePts),
                                stringsAsFactors = F)
  p <- ggplot2::ggplot(head(dat.plot.eigenv,n=30),mapping = aes(PC,eigenv)) +
    geom_point(aes(colour=isKneePts),show.legend=F) +
    scale_colour_manual(values = c("TRUE"="#E41A1C","FALSE"="#377EB8")) +
    theme_bw()
  print(p)
}

