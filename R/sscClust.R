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
  }
  if(is.null(rowData(obj)[["display.name"]])){
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
#' @importFrom scran trendVar decomposeVar
#' @param obj object of SingleCellExperiment
#' @param method method to be used, can be one of "sd", mean.sd, trendVar. (default: "sd")
#' @param sd.n top number of genes (default 1500)
#' @param mean.thre numeric; threshold for mean, used in trendVar method (default 0.1)
#' @param assay.name which assay to be used (default "exprs")
#' @param var.block character; specify the uninteresting factors by formula. E.g. "~patient" (default NULL)
#' @param reuse logical; don't calculate if the query is already available. (default: F)
#' @param out.prefix character; if not NULL, output prefix. (default: F)
#' @details Method "sd": calculate the standard deviation of each gene and sort decreasingly, the top `sd.n` genes are the
#' variable genes. Method "trendVar": fit the trend between variance and mean, and decompose each gene's variance into
#' 'tech' part(fitted value) and 'bio' part (residual value), then select genes according FDR and mean threshold. Note,
#' when using "trendVar", will use expression data stored in "norm_exprs" slot of `obj`, no matter what `assay.name` is.
#' @return an object of \code{SingleCellExperiment} class
#' @export
ssc.variableGene <- function(obj,method="sd",sd.n=1500,mean.thre=0.1,assay.name="exprs",
                             var.block=NULL,
                             reuse=F,out.prefix=NULL)
{
  if(method=="sd")
  {
    if(!reuse || is.null(metadata(obj)$ssc[["variable.gene"]][["sd"]])){
      row.sd <- apply(assay(obj,assay.name),1,sd)
      metadata(obj)$ssc[["variable.gene"]][["sd"]] <- names(head(row.sd[order(row.sd,decreasing = T)],n=sd.n))
    }
  }else if(method=="mean.sd")
  {
    if(!reuse || is.null(metadata(obj)$ssc[["variable.gene"]][["mean.sd"]])){
      row.sd <- apply(assay(obj,assay.name),1,sd)
      row.mean <- apply(assay(obj,assay.name),1,mean)
      info.gene.sd <- row.sd[row.sd>1 & row.mean>1]
      info.gene.sd <- sort(info.gene.sd,decreasing = T)
      metadata(obj)$ssc[["variable.gene"]][["mean.sd"]] <- names(head(info.gene.sd,n=sd.n))
    }
  }else if(method=="trendVar")
  {
    trendVar.min.mean <- 0.1
    if(!is.null(var.block)){
      var.design <- model.matrix(as.formula(var.block),data = colData(obj))
    }else{
      var.design <- NULL
    }
    var.fit <- scran::trendVar(obj, parametric=TRUE, span=0.3, min.mean=trendVar.min.mean,
                               design=var.design,
                               use.spikes=F, assay.type="norm_exprs")
                               ####use.spikes=F, assay.type=assay.name)
    var.out <- scran::decomposeVar(obj, var.fit, assay.type="norm_exprs")
    ####var.out <- scran::decomposeVar(obj, var.fit, assay.type=assay.name)
    var.out <- var.out[order(var.out$FDR,-var.out$bio/var.out$total),]
    f.var <- var.out$FDR<0.001 & var.out$mean>mean.thre
    metadata(obj)$ssc[["variable.gene"]][["trendVar"]] <- rownames(var.out)[f.var]
    #head(var.out)
    ### debug
    row.sd <- apply(assay(obj,assay.name),1,sd)
    sd.thre <- sort(row.sd,decreasing = T)[min(sd.n,length(row.sd))]
    print("debug info (ssc.variableGene)")
    list.a <- names(row.sd)[row.sd>=sd.thre]
    list.b <- rownames(var.out)[f.var]
    print(length(list.a))
    print(length(list.b))
    print(length(intersect(list.a,list.b)))
    if(!is.null(out.prefix) && file.exists(dirname(out.prefix))){
      pdf(sprintf("%s.varGene.%s.pdf",out.prefix,method),width = 10,height = 4)
      opar=par(mfcol=c(1,2))
      plot(var.out$mean, var.out$total,pch=16,cex=0.3, xlab="Mean log-expression",
           ylab="Variance of log-expression")
      points(var.out$mean[f.var],
             var.out$total[f.var], col="red", pch=16,cex=0.3)
      curve(var.fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
      plot(var.out$mean, var.out$total,pch=16,cex=0.3, xlab="Mean log-expression",
           ylab="Variance of log-expression")
      points(var.out[names(row.sd)[row.sd>=sd.thre],"mean"],
             var.out[names(row.sd)[row.sd>=sd.thre],"total"], col="red", pch=16,cex=0.3)
      curve(var.fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
      par(opar)
      dev.off()
    }
    ### end of debug
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
#' @param method character; method to be used for dimension reduction, should be one of (pca, tsne, iCor, zinbwave). (default: "iCor")
#' @param method.vgene character; method to identify variable genes. (default: sd)
#' @param pca.npc integer; number of pc be used. Only for reduction method "pca". (default: NULL)
#' @param tSNE.usePCA logical; whether use PCA before tSNE. Only for reduction method "tsne". (default: T)
#' @param tSNE.perplexity logical; perplexity parameter. Used in all Rtsne() calling. (default: 30)
#' @param autoTSNE logical; Wheter generate automatically a tSNE map when reduction method is "pca" or "iCor". (default: T)
#' @param dim.name character; store the reduced data under the name in the obj's reducedDim SimpleList. If i
#' it is NULL, infer from method. (default: NULL)
#' @param iCor.niter integer; number of iteration of calculating the correlation. Used in reduction method "iCor". (default: 1)
#' @param iCor.method character; method to calculate correlation between samples,
#' should be one of "spearman" and "pearson". (default "spearman")
#' @param zinbwave.K integer, zinbwave parameter, number of latent variables. (default: 20)
#' @param zinbwave.X character, zinbwave parameter, cell-level covariates. (default: "~patient")
#' @param reuse logical; don't calculate if the query is already available. (default: F)
#' @param seed integer; seed of random number generation. (default: NULL)
#' @param ncore integer; nuber of CPU cores to use. if NULL, automatically detect the number. (default: NULL)
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
                          zinbwave.K=20, zinbwave.X="~patient",
                          autoTSNE=T,
                          dim.name=NULL,
                          iCor.niter=1,iCor.method="spearman",
                          reuse=F,ncore=NULL,seed=NULL)
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
  if(is.null(dim.name)){ dim.name <- method }
  vgene <- metadata(obj)$ssc[["variable.gene"]][[method.vgene]]
  vgene <- intersect(vgene,rownames(obj))
  if(!reuse || !(dim.name %in% reducedDimNames(obj)) ){
    if(is.null(seed)){
      myseed <- as.integer(Sys.time())
    }else{
      myseed <- seed
    }
    loginfo(sprintf("set.seed(%s) for ssc.reduceDim\n",as.character(myseed)))
    set.seed(myseed)

    if(method=="pca"){
        pca.res <- prcomp(t(assay(obj[vgene,],assay.name)))
        ### find elbow point and get number of components to be used
        pca.res$eigenv.prop <- pca.res$sdev/sum(pca.res$sdev)
        pca.res$eigengap <- sapply(seq_len(length(pca.res$eigenv.prop)-1),function(i){ pca.res$eigenv.prop[i]-pca.res$eigenv.prop[i+1] })
        ### method 1
        ######pca.res$kneePts <- which(pca.res$eigengap<1e-4)[1]
        ### method 2 (max distance to the line defined by the first and last point in the scree plot)
        pca.res$kneePts <- findKneePoint(head(pca.res$eigenv.prop,n=100))
        if(is.null(pca.npc)){
          if(!is.na(pca.res$kneePts)){ pca.npc <- pca.res$kneePts
          }else{ pca.npc <- 30 }
        }
        pca.npc <- min(pca.npc,ncol(pca.res$x))
        pca.res$npc <- pca.npc
        loginfo(sprintf("set pca.npc to %d while kneePts is at %d (ssc.reduceDim)\n",pca.npc,
                        if(!is.na(pca.res$kneePts)) pca.res$kneePts else -1))
        ### save to object
        metadata(obj)$ssc$pca.res <- pca.res
        proj_data <- pca.res$x[,1:pca.npc,drop=F]
        if(autoTSNE){ reducedDim(obj,sprintf("%s.tsne",dim.name)) <- run.tSNE(proj_data,tSNE.usePCA=F,tSNE.perplexity) }
    }else if(method=="tsne"){
      proj_data <- run.tSNE(t(assay(obj[vgene,],assay.name)),tSNE.usePCA,tSNE.perplexity)
    }else if(method=="iCor"){
      proj_data <- assay(obj[vgene,],assay.name)
      while(iCor.niter>0){
        ##proj_data <- cor(proj_data,method=iCor.method)
        proj_data <- cor.BLAS(t(proj_data),method=iCor.method,nthreads = ncore)
        iCor.niter <- iCor.niter-1
      }
      if(autoTSNE) { reducedDim(obj,sprintf("%s.tsne",dim.name)) <- run.tSNE(proj_data,tSNE.usePCA=F,tSNE.perplexity) }
    }else if(method=="zinbwave"){
      res.zinb <- run.zinbWave(obj,assay.name=assay.name,vgene=vgene,n.cores=ncore,
                               zinbwave.K=zinbwave.K, zinbwave.X=zinbwave.X,verbose=F)
      proj_data <- getW(res.zinb)
      colnames(proj_data) <- sprintf("W%d",seq_len(ncol(proj_data)))
      if(autoTSNE) { reducedDim(obj,sprintf("%s.tsne",dim.name)) <- run.tSNE(proj_data,tSNE.usePCA=F,tSNE.perplexity) }
    }
    reducedDim(obj,dim.name) <- proj_data
  }
  return(obj)
}


#' Perform clustering using reduced data
#' @importFrom factoextra eclust
#' @importFrom stats kmeans dist
#' @importFrom ADPclust adpclust
#' @importFrom densityClust densityClust findClusters
#' @importFrom scran buildSNNGraph
#' @importFrom igraph cluster_fast_greedy cluster_leading_eigen cluster_infomap
#' cluster_label_prop cluster_louvain cluster_optimal cluster_spinglass cluster_walktrap
#' cluster_edge_betweenness
#' @importFrom plyr llply
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param method.reduction character; which dimention reduction method to be used, should be one of
#' "iCor", "pca" and "none". (default: "iCor")
#' @param method character; clustering method to be used, should be one of "kmeans", "hclust", "SNN", "dpclust", "adpclust" and "SC3". (default: "kmeans")
#' @param k.batch integer; number of clusters to be evaluated. (default: 2:6)
#' @param method.vgene character; variable gene identification method used. (default: "sd")
#' @param SNN.k integer; number of shared NN. (default: 10)
#' @param SNN.method character; cluster method applied on SNN， one of "greedy", "eigen", "infomap",
#' "prop", "louvain", "optimal", "spinglass", "walktrap", "betweenness". (default: "eigen")
#' @param SC3.biology logical, SC3 parameter, whether calcualte biology. (default: T)
#' @param SC3.markerplot.width integer, SC3 parameter, with of the marker plot (default: 15)
#' @param dpclust.rho numberic; cuttoff of rho, if it is NULL, infer frome the data (default: NULL)
#' @param dpclust.delta numberic; cuttoff of delta, if it is NULL, infer frome the data (default: NULL)
#' @param out.prefix character; output prefix, if not NULL, some plots of intermediate result will be produced. (default: NULL)
#' @param parlist list; if not NULL, use th parameters in it. (default: NULL)
#' @param seed integer; seed of random number generation. (default: NULL)
#' @param ncore integer; nuber of CPU cores to use. if NULL, automatically detect the number. (default: NULL)
#' @details If no dimension reduction performed or method is "none", expression data of variable genes,
#' which can be speficed by method.vgene, will be used for clustering. Otherwise, the reduced data specified by
#' method.reduction will be used. The cluster label will stored in the colData of the object
#' of \code{singleCellExperiment} class, with colname in the format of \{method.reduction\}.\{method\}k\{k\}
#' where \{k\} get value(s) from k.batch.
#' @return an object of \code{SingleCellExperiment} class with cluster labels added.
#' @export
ssc.clust <- function(obj, assay.name="exprs", method.reduction="iCor",
                      method="kmeans", k.batch=2:6,
                      method.vgene="sd",
                      SNN.k=10,SNN.method="eigen",
                      SC3.biology=T,SC3.markerplot.width=15,
                      dpclust.rho=NULL,dpclust.delta=NULL,
                      parlist=NULL,
                      out.prefix=NULL,seed=NULL,ncore=NULL)
{
  clust.res <- NULL
  res.list <- list()
  ### check transformed data
  if(method.reduction=="none"){
    warning(sprintf("The dimention reduction should be performed!"))
    vgene <- metadata(obj)$ssc[["variable.gene"]][[method.vgene]]
    dat.transformed <- t(assay(obj[vgene,],assay.name))
  }else if( (!method.reduction %in% reducedDimNames(obj)) ){
    dat.transformed <- NULL
  }else{
    dat.transformed <- reducedDim(obj,method.reduction)
  }
  if(is.null(dat.transformed) && method!="SC3"){
    warning("dat.transformed is null !!")
    metadata(obj)$ssc$clust.res[[method]] <- NULL
    if(method %in% c("adpclust","dpclust","SNN")){
      k <- "auto"
      colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%d",rep(0,ncol(obj)))
    }else{
      for(k in k.batch){
        colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%d",rep(0,ncol(obj)))
      }
    }
    ##print(obj)
    return(obj)
  }
  ### check method

  ###
  if(is.null(seed)){
    myseed <- as.integer(Sys.time())
  }else{
    myseed <- seed
  }
  loginfo(sprintf("set.seed(%s) for ssc.clust\n",as.character(myseed)))
  set.seed(myseed)

  ### No k needed for methods: "adpclust","dpclust","SNN"
  if(method=="adpclust"){
    clust.res <- ADPclust::adpclust(dat.transformed,nclust = k.batch)
    k <- "auto"
    colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%d",clust.res$clusters)
    if(!is.null(out.prefix)){
      dir.create(dirname(out.prefix),showWarnings = F,recursive = T)
      pdf(sprintf("%s.adpclust.diagnostic.pdf",out.prefix),width=11,height = 5)
      plot(clust.res)
      dev.off()
      ssc.plot.tsne(obj,plotDensity = T,reduced.name = sprintf("%s",method.reduction),
                    peaks = clust.res$centers,
                    out.prefix = sprintf("%s.aadpclust",out.prefix),base_aspect_ratio = 1.4)
    }
    res.list[[as.character(k)]] <- clust.res
  }else if(method=="dpclust"){
    dist.obj <- stats::dist(dat.transformed)
    clust.res <- densityClust::densityClust(dist.obj, gaussian = T)
    if(is.null(dpclust.rho)){ dpclust.rho <- quantile(clust.res$rho, probs = 0.90) }
    if(is.null(dpclust.delta)){ dpclust.delta <- quantile(clust.res$delta, probs = 0.95) }
    ## overwritet the parameter using those in parlist
    if(!is.null(parlist)){
      if("rho" %in% names(parlist)){ dpclust.rho <- parlist[["rho"]] }
      if("delta" %in% names(parlist)){ dpclust.delta <- parlist[["delta"]] }
    }
    cat(sprintf("(dpclust.rho: %4.2f)\n",dpclust.rho))
    cat(sprintf("(dpclust.delta: %4.2f)\n",dpclust.delta))
    if(sum(clust.res$rho >= dpclust.rho & clust.res$delta >= dpclust.delta)<2){
      clust.res$clusters <- rep(1,length(clust.res$rho))
      clust.res$peaks <- NA
    }else{
      clust.res <- densityClust::findClusters(clust.res,rho = dpclust.rho,
                                              delta = dpclust.delta,plot=F)
    }
    k <- "auto"
    colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%d",clust.res$clusters)
    res.list[[as.character(k)]] <- clust.res
    if(!is.null(out.prefix)){
      pdf(sprintf("%s.decision.pdf",out.prefix),width = 5,height = 5)
      plot(clust.res$rho, clust.res$delta, main = 'Decision graph', xlab = expression(rho),
           ylab = expression(delta))
      if (all(!is.na(clust.res$peaks))) {
        points(clust.res$rho[clust.res$peaks], clust.res$delta[clust.res$peaks],
               col = 2:(1 + length(clust.res$peaks)),pch = 19)
      }
      abline(v=dpclust.rho,lty=2)
      abline(h=dpclust.delta,lty=2)
      plot(sort(clust.res$rho * clust.res$delta,decreasing = T),ylab="rho * delta")
      dev.off()
      ssc.plot.tsne(obj,plotDensity = T,reduced.name = sprintf("%s",method.reduction),
                    peaks = if(all(!is.na(clust.res$peaks))) clust.res$peaks else NULL,
                    out.prefix = sprintf("%s.dpclust",out.prefix),base_aspect_ratio = 1.4)
    }
  }else if(method=="SNN"){
    if(!is.null(metadata(obj)$ssc$clust.res[["snn.gr"]])){
      snn.gr <- metadata(obj)$ssc$clust.res[["snn.gr"]]
    }else{
      loginfo(sprintf("buildSNNGraph with k=%d\n",SNN.k))
      snn.gr <- scran::buildSNNGraph(t(dat.transformed), k=SNN.k,d=NA)
    }
    loginfo(sprintf("cluster using method %s begin\n",SNN.method))
    if(SNN.method=="greedy"){
      clust.res <- igraph::cluster_fast_greedy(snn.gr)
    }else if(SNN.method=="eigen"){
      clust.res <- igraph::cluster_leading_eigen(snn.gr)
      ##clust.res <- igraph::cluster_leading_eigen(snn.gr,weights = igraph::E(snn.gr)$weight)
    }else if(SNN.method=="infomap"){
      clust.res <- igraph::cluster_infomap(snn.gr)
    }else if(SNN.method=="prop"){
      clust.res <- igraph::cluster_label_prop(snn.gr)
    }else if(SNN.method=="louvain"){
      clust.res <- igraph::cluster_louvain(snn.gr)
    }else if(SNN.method=="optimal"){
      clust.res <- igraph::cluster_optimal(snn.gr)
    }else if(SNN.method=="spinglass"){
      clust.res <- igraph::cluster_spinglass(snn.gr)
    }else if(SNN.method=="walktrap"){
      clust.res <- igraph::cluster_walktrap(snn.gr)
    }else if(SNN.method=="betweenness"){
      clust.res <- igraph::cluster_edge_betweenness(snn.gr)
    }
    metadata(obj)$ssc$clust.res[["snn.gr"]] <- snn.gr
    loginfo(sprintf("cluster using method %s finished\n",SNN.method))
    k <- "auto"
    colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%d",clust.res$membership)
    res.list[[as.character(k)]] <- clust.res
  }else if(method=="SC3"){
    vgene <- metadata(obj)$ssc[["variable.gene"]][[method.vgene]]
    obj.tmp <- run.SC3(obj[vgene,],assay.name = assay.name,out.prefix=out.prefix,n.cores = ncore,ks=k.batch,
                       SC3.biology=SC3.biology,SC3.markerplot.width=SC3.markerplot.width)
    colData(obj) <- colData(obj.tmp)
    .cls.labbel <- grep("sc3_\\d+_clusters",names(colData(obj)),perl=T,value = T)
    for(.cls.l in .cls.labbel){
      colData(obj)[[.cls.l]] <- as.character(colData(obj)[[.cls.l]])
    }
    metadata(obj)$sc3 <- metadata(obj.tmp)$sc3
    for(.dname in names(metadata(obj)$sc3$transformations)){
      reducedDim(obj,sprintf("sc3.%s",.dname)) <- metadata(obj)$sc3$transformations[[.dname]]
      reducedDim(obj,sprintf("sc3.%s.tsne",.dname)) <- run.tSNE(reducedDim(obj,sprintf("sc3.%s",.dname)),tSNE.usePCA=F,30)
    }
    res.list[["sc3.biology"]] <- rowData(obj.tmp)[,grepl("^(sc3_|display.name|feature_symbol)",names(rowData(obj.tmp)),perl = T),drop=F]
  }else{
    ### for other methods need k, llply on k.batch
    RhpcBLASctl::omp_set_num_threads(1)
    registerDoParallel(cores = ncore)
    .tmp.ret <- NULL
    tryCatch({
      .tmp.ret <- llply(k.batch,function(k){
        if(method=="kmeans"){
            loginfo(sprintf("begin kmeans(..,k=%d) ...",k))
            clust.res <- kmeans(dat.transformed,k,iter.max=1000,nstart=50)
            loginfo(sprintf("end kmeans(..,k=%d) ...",k))
        }else if(method=="hclust"){
            clust.res <- factoextra::eclust(dat.transformed, "hclust", k = k, method = "complete", graph = FALSE)
        }
        list("clust.res"=clust.res,"label"=sprintf("C%d",clust.res$cluster))
        #colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%d",clust.res$cluster)
      },.progress = "none",.parallel=T)
    },error=function(e){
      cat(sprintf("Error occur in ssc.clust --> llply(k.batch,...).\n"))
      print(e)
    })
    names(.tmp.ret) <- sprintf("%d",k.batch)
    for(k in k.batch){
        colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- .tmp.ret[[as.character(k)]][["label"]]
        res.list[[as.character(k)]] <- .tmp.ret[[as.character(k)]][["clust.res"]]
    }
  }
  metadata(obj)$ssc$clust.res[[method]] <- res.list
  return(obj)
}

#' Clustering with subsampling and classification
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom stats cor
#' @importFrom plyr llply
#' @importFrom MASS ginv
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
#' @param vis.proj logical; whether get low dimensional representation for visualization. (default: F)
#' @param ncore integer; number of cpu to use. (default: NULL)
#' @param k.batch integer; number of clusters to be evaluated. (default: 2:6)
#' @param seed integer; seed of random number generation. (default: NULL)
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
                                               vis.proj=F,
                                               ncore=NULL,
                                               k.batch=2:6,seed=NULL)
{
  #### subsampling
  if(frac>1){
      n.sub <- floor(frac)
  }else{
      n.sub <- floor(ncol(obj)*frac)
  }
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
  metadata(obj.train)$ssc[["variable.gene"]][[method.vgene]] <- intersect(metadata(obj.train)$ssc[["variable.gene"]][[method.vgene]],
                                                                        rownames(obj.train))
  metadata(obj.pred)$ssc[["variable.gene"]][[method.vgene]] <- intersect(metadata(obj.pred)$ssc[["variable.gene"]][[method.vgene]],
                                                                         rownames(obj.pred))

  reducedDims(obj.train) <- SimpleList()
  #### reduction
  if(method.reduction=="none"){
    warning(sprintf("Reducing the dimensions first is recommended"))
  }else if(method.reduction %in% c("pca","iCor")){
    loginfo(sprintf("begin obj.train=ssc.reduceDim(,...,method=%s) ...",method.reduction))
    obj.train <- ssc.reduceDim(obj.train,assay.name = assay.name,method=method.reduction,
                              method.vgene=method.vgene,
                              pca.npc=pca.npc,autoTSNE = F,seed = seed, ### disable auto tSNE temporarily
                              iCor.niter=iCor.niter,iCor.method="spearman",ncore=ncore)
    loginfo(sprintf("end obj.train=ssc.reduceDim(,...,method=%s) ...",method.reduction))
    if(!all(rownames(obj.train)==rownames(obj.pred))){
      stop("The genes of obj.train and obj.pred are different!")
    }
    vgene <- metadata(obj.pred)$ssc$variable.gene[[method.vgene]]

    if(method.reduction=="iCor"){
      loginfo(sprintf("begin A=cor.BLAS() ..."))
      A <- cor.BLAS(t(assay(obj.pred[vgene,],assay.name)),t(assay(obj.train[vgene,],assay.name)),method = "spearman",nthreads=ncore)
      loginfo(sprintf("end A=cor.BLAS() ..."))
      #### for (tsne) visualization
      ###
      #vis.proj <- F
      if(vis.proj){
          S <- cor.BLAS(t(assay(obj.train[vgene,],assay.name)),method = "spearman",nthreads=ncore)
          dat.map.train <- reducedDim(obj.train,sprintf("%s.tsne",method.reduction))
          ### S may be computationally singular
          ####dat.map.pred <- A %*% solve(S) %*% dat.map.train
          dat.map.pred <- A %*% MASS::ginv(S) %*% dat.map.train
          reducedDim(obj.pred,sprintf("%s.tsne",method.reduction)) <- dat.map.pred
          dat.map.obj <- rbind(dat.map.train,dat.map.pred)
          reducedDim(obj,sprintf("%s.tsne",method.reduction)) <- dat.map.obj[colnames(obj),,drop=F]
#      }else{
#          obj <- ssc.reduceDim(obj,assay.name = assay.name,method=method.reduction,
#                              method.vgene=method.vgene,
#                              pca.npc=pca.npc,autoTSNE = T,seed = seed,
#                              iCor.niter=iCor.niter,iCor.method="spearman")
      }
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
  obj.train <- ssc.clust(obj.train,assay.name = assay.name, method.reduction=method.reduction,
                         method=method.clust, k.batch=k.batch,seed = seed,ncore=ncore)
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
  RhpcBLASctl::omp_set_num_threads(1)
  registerDoParallel(cores = ncore)
  .tmp.ret <- NULL
  tryCatch({
    .tmp.ret <- llply(k.batch,function(kk){
        cls.label <- colData(obj.train)[,sprintf("%s.%s.k%s",method.reduction,method.clust,kk)]
        names(cls.label) <- colnames(obj.train)
        loginfo(sprintf("begin run.%s(..) base on clustering result using k=%d ...",toupper(method.classify),kk))
        if(method.classify=="RF") {
          res.pred <- run.RF(data.train, as.factor(cls.label), data.pred,do.norm = T)
        }else if(method.classify=="knn") {
          res.pred <- run.KNN(data.train,as.factor(cls.label), data.pred,k=1)
        }else if(method.classify=="svm") {
          res.pred <- run.SVM(data.train,as.factor(cls.label),data.pred,kern="linear")
        }
        loginfo(sprintf("end run.%s(..) base on clustering result using k=%d ...",toupper(method.classify),kk))
        pred.label <- structure(as.character(res.pred$ylabel), names=names(res.pred$ylabel))
        ret.label <- c(cls.label,pred.label)
        list("ret.label"=ret.label)
    },.progress = "none",.parallel=T)
  },error=function(e){
      cat(sprintf("Error occur in ssc.clust --> llply(k.batch,...).\n"))
      print(e)
  })
  names(.tmp.ret) <- sprintf("%d",k.batch)
  for(kk in k.batch){
    ret.label <- .tmp.ret[[as.character(kk)]][["ret.label"]]
    colData(obj)[names(ret.label),sprintf("%s.%s.k%s",method.reduction,method.clust,kk)] <- ret.label
  }
  return(obj)
}

#' Wrapper for running all the pipeline
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param method.vgene character; variable gene identification method used. (default: "sd")
#' @param mean.thre numeric; threshold for mean, used in trendVar method (default 0.1)
#' @param var.block character; specify the uninteresting factors by formula. E.g. "~patient".
#' used in trendVar method (default NULL)
#' @param sd.n integer; top number of genes as variable genes (default 1500)
#' @param de.n integer; number of differential genes used for refined geneset for another run of clustering (default 1500)
#' @param method.reduction character; which dimention reduction method to be used, should be one of
#' "iCor", "pca", "zinbwave" and "none". (default: "iCor")
#' @param method.clust character; clustering method to be used, should be one of "kmeans", "hclust", "SNN", "adpclust" and "SC3. (default: "kmeans")
#' @param method.classify character; method used for classification, one of "knn" and "RF". (default: "knn")
#' @param pca.npc integer; number of pc be used. Only for reduction method "pca". (default: NULL)
#' @param iCor.niter integer; number of iteration of calculating the correlation. Used in reduction method "iCor". (default: 1)
#' @param iCor.method character; correlation method, one of "spearman", "pearson" (default: "spearman")
#' @param zinbwave.K integer, zinbwave parameter, number of latent variables. (default: 20)
#' @param zinbwave.X character, zinbwave parameter, cell-level covariates. (default: "~patient")
#' @param subsampling logical; whether cluster using the subsampling->cluster->classification method. (default: F)
#' @param sub.frac numeric; subsample to frac of original samples. (default: 0.4)
#' @param sub.use.proj logical; whether use the projected data for classification. (default: T)
#' @param sub.vis.proj logical; whether get low dimensional representation for visualization, only used in downsample mode. (default: F)
#' @param k.batch integer; number of clusters to be evaluated. (default: 2:6)
#' @param refineGene logical; whether perform second round demension reduction and clustering pipeline using the differential
#' genes found by the first round cluster result. (default: F)
#' @param nIter integer; number of iterative clustering in sub-cluster. (default: 1)
#' @param do.DE logical; perform DE analysis when clustering finished. (default: F)
#' @param out.prefix character; output prefix, if not NULL, some plots of intermediate result will be produced. (default: NULL)
#' @param parfile character; parameter files, if not NULL, will use the settings. must contain a list named
#' `parlist`. (default: NULL)
#' @param ncore integer; nuber of CPU cores to use. if NULL, automatically detect the number. (default: NULL)
#' @param reuse logical; don't calculate if the query is already available. (default: F)
#' @param seed integer; seed of random number generation. (default: NULL)
#' @param ... parameters pass to clustering methods
#' @details run the pipeline of variable gene identification, dimension reduction, clustering.
#' @seealso \code{\link{ssc.variableGene}} for variable genes' identification, \code{\link{ssc.reduceDim}}
#' for dimension reduction, \code{\link{ssc.clust}} for clustering using all data
#' and \code{\link{ssc.clustSubsamplingClassification}} for clustering with subsampling.
#' @return an object of \code{SingleCellExperiment} class with cluster labels added.
#' @export
ssc.run <- function(obj, assay.name="exprs",
                    method.vgene="sd",
                    sd.n=1500,
                    mean.thre=0.1,
                    var.block=NULL,
                    method.reduction="iCor",
                    method.clust="kmeans",
                    method.classify="knn",
                    pca.npc=NULL,
                    iCor.niter=1,
                    iCor.method="spearman",
                    zinbwave.K=20, zinbwave.X="~patient",
                    subsampling=F,
                    sub.frac=0.4,
                    sub.use.proj=T,
                    sub.vis.proj=F,
                    k.batch=2:6,
                    refineGene=F,
                    de.n=1500,
                    nIter=1,
                    out.prefix=NULL,
                    parfile=NULL,
                    reuse=F,
                    ncore=NULL,
                    seed=NULL,
                    do.DE=F,...)
{
  ### some checking
  if(is.null(colnames(obj))){
    stop("colnames of obj is NULL!!!")
  }
  if(is.null(rownames(obj))){
    stop("rownames of obj is NULL!!!")
  }
  #rid <- "11"
  #obj <- ssc.variableGene(obj,method=method.vgene,sd.n=sd.n,assay.name=assay.name)
  if(!subsampling){
    if(!is.null(parfile) && file.exists(parfile)){
      source(parfile)
      cat(sprintf("parfile: %s\n",parfile))
      print(parlist)
    }else{
      parlist <- NULL
    }
    runOneIter <- function(obj,rid,k.batch,level=1){
      if(!is.null(parlist) && rid %in% names(parlist)){
        parlist.rid <- parlist[[rid]]
        if("k" %in% names(parlist.rid) && parlist.rid[["k"]]==1){
          return(obj)
        }
      }else{
        parlist.rid <- NULL
      }
      loginfo(sprintf("select variable genes ... (%s)",rid))
      obj <- ssc.variableGene(obj,method=method.vgene,sd.n=sd.n,mean.thre = mean.thre,
                              assay.name=assay.name,reuse = reuse,var.block = var.block,
                              out.prefix = sprintf("%s.%s",out.prefix,rid))
      loginfo(sprintf("reduce dimensions ... (%s)",rid))
      obj <- ssc.reduceDim(obj,assay.name=assay.name,
                           method=method.reduction,
                           pca.npc = pca.npc,
                           iCor.niter = iCor.niter,
                           iCor.method = iCor.method,
                           zinbwave.K = zinbwave.K, zinbwave.X = zinbwave.X,
                           method.vgene=method.vgene,
                           ncore = ncore,
                           seed = seed,
                           reuse = reuse)
      loginfo(sprintf("clustering ... (%s)",rid))
      obj <- ssc.clust(obj, assay.name=assay.name,
                       method.reduction=if(method.clust %in% c("adpclust","dpclust")) sprintf("%s.tsne",method.reduction) else method.reduction,
                       method=method.clust, k.batch=k.batch,
                       out.prefix = if(is.null(out.prefix)) NULL else sprintf("%s.%s",out.prefix,rid),
                       seed = seed,
                       method.vgene=method.vgene, parlist = parlist.rid, ...)

      .xlabel <- NULL
      if(method.clust %in% c("adpclust","dpclust")){
        .xlabel <- sprintf("%s.tsne.%s.kauto",method.reduction,method.clust)
      }else if(method.clust=="SNN"){
        .xlabel <- sprintf("%s.%s.kauto",method.reduction,method.clust)
      }
      if(!is.null(out.prefix) && !is.null(.xlabel)){
        ssc.plot.tsne(obj,columns = c(.xlabel),
                      reduced.name = if(method.clust %in% c("adpclust","dpclust")) sprintf("%s.tsne",method.reduction) else method.reduction,
                      out.prefix = sprintf("%s.%s",out.prefix,rid),
                      base_aspect_ratio = 1.4)
      }
      ### other method need determine the best k. not implemented yet.
      if(refineGene && method.clust %in% c("adpclust","dpclust","SNN")){
        do.secondRun <- T
        if(!is.null(parlist) && sprintf("%s.de",rid) %in% names(parlist)){
          parlist.rid <- parlist[[sprintf("%s.de",rid)]]
          if("k" %in% names(parlist.rid) && parlist.rid[["k"]]==1){
            do.secondRun <- F
          }
        }else{
          parlist.rid <- NULL
        }
        if(do.secondRun)
        {
          ### adpclust automatically use tsne data
          loginfo(sprintf("find DE genes ... (%s)",rid))
          de.out <- findDEGenesByAOV(xdata = assay(obj,assay.name),
                                     xlabel = colData(obj)[,.xlabel],
                                     n.cores = ncore,
                                     gid.mapping = rowData(obj)[,"display.name"])
          if(!is.null(de.out) && nrow(de.out$aov.out.sig)>30){
            metadata(obj)$ssc[["de.res"]][[rid]] <- de.out
            metadata(obj)$ssc[["variable.gene"]][["refine.de"]] <- head(de.out$aov.out.sig$geneID,n=de.n)
            loginfo(sprintf("reduce dimensions using DE genes ... (%s)",rid))
            obj <- ssc.reduceDim(obj,assay.name=assay.name,
                         method=method.reduction,
                         pca.npc = pca.npc,
                         iCor.niter = iCor.niter,
                         iCor.method = iCor.method,
                         zinbwave.K = zinbwave.K, zinbwave.X = zinbwave.X,
                         ncore = ncore,
                         seed = seed,
                         dim.name = sprintf("de.%s",method.reduction),
                         method.vgene="refine.de",reuse = reuse)
            loginfo(sprintf("clust using DE genes ... (%s)",rid))
            obj <- ssc.clust(obj, assay.name=assay.name,
                             method.reduction=if(method.clust %in% c("adpclust","dpclust")) sprintf("de.%s.tsne",method.reduction) else sprintf("de.%s",method.reduction),
                             ##method.reduction=if(method.clust %in% c("adpclust","dpclust")) sprintf("%s.tsne",method.reduction) else method.reduction,
                             method=method.clust, k.batch=k.batch,
                             seed = seed,
                             out.prefix = if(is.null(out.prefix)) NULL else sprintf("%s.%s.refineG",out.prefix,rid),
                             method.vgene="refine.de", parlist = parlist.rid, ...)

            if(method.clust %in% c("adpclust","dpclust")){
              colData(obj)[,.xlabel] <- colData(obj)[,sprintf("de.%s.tsne.%s.kauto",method.reduction,method.clust)]
              colData(obj)[,sprintf("de.%s.tsne.%s.kauto",method.reduction,method.clust)] <- NULL
            }else if(method.clust=="SNN"){
              colData(obj)[,.xlabel] <- colData(obj)[,sprintf("de.%s.%s.kauto",method.reduction,method.clust)]
              colData(obj)[,sprintf("de.%s.%s.kauto",method.reduction,method.clust)] <- NULL
            }
            if(!is.null(out.prefix) && !is.null(.xlabel)){
              ssc.plot.tsne(obj,columns = c(.xlabel),
                            reduced.name = if(method.clust %in% c("adpclust","dpclust")) sprintf("de.%s.tsne",method.reduction) else sprintf("de.%s",method.reduction),
                            out.prefix = sprintf("%s.%s.refineG",out.prefix,rid),
                            base_aspect_ratio = 1.4)
            }
          }else{
            warning("The number of DE genes is less than 30, NO second round clustering using DE genes will be performed!!")
          }
        }
      }
      if(level<nIter){
        clustLabel <- if(method.clust %in% c("adpclust","dpclust")) sprintf("%s.tsne.%s.kauto",method.reduction,method.clust) else sprintf("%s.%s.kauto",method.reduction,method.clust)
        if(clustLabel %in% colnames(colData(obj))){
          clusterNames <- sort(unique(colData(obj)[,clustLabel]))
          for(cls in clusterNames){
            f.cls <- colData(obj)[,clustLabel]==cls
            obj.cls <- obj[,f.cls]
            cls.rid <- sprintf("%s.L%s%s",rid,level+1,cls)
            ### clear the metadata of obj.cls
            obj.cls <- ssc.build(obj.cls)
            if(ncol(obj.cls)>10)
            {
              nsamples <- ncol(obj.cls)
              obj.cls <- runOneIter(obj.cls,cls.rid,k.batch = k.batch[k.batch < nsamples],level+1)
              # update cluster's label
              colData(obj)[f.cls,clustLabel] <- sprintf("%s.L%s%s",cls.rid,level+2,
                                                        colData(obj.cls)[,clustLabel])
            }else{
              # update cluster's label
              colData(obj)[f.cls,clustLabel] <- sprintf("%s.L%s%s",cls.rid,level+2,"C1")
            }
          }
        }
      }
      return(obj)
    }
    obj <- runOneIter(obj,"L1C1",k.batch = k.batch)
    if(do.DE && method.clust %in% c("adpclust","dpclust","SNN")){
      .xlabel <- NULL
      if(method.clust %in% c("adpclust","dpclust")){
        .xlabel <- colData(obj)[,sprintf("%s.tsne.%s.kauto",method.reduction,method.clust)]
      }else if(method.clust=="SNN"){
        .xlabel <- colData(obj)[,sprintf("%s.%s.kauto",method.reduction,method.clust)]
      }
      de.out <- findDEGenesByAOV(xdata = assay(obj,assay.name),
                                 xlabel = .xlabel,
                                 gid.mapping = rowData(obj)[,"display.name"])
      metadata(obj)$ssc[["de.res"]][["L1C1"]] <- de.out
      metadata(obj)$ssc[["variable.gene"]][["de"]] <- head(de.out$aov.out.sig$geneID,n=sd.n)
      ### for general visualization
      obj <- ssc.reduceDim(obj,assay.name=assay.name, method="tsne",
                           zinbwave.K = zinbwave.K, zinbwave.X = zinbwave.X,
                           method.vgene="de",dim.name = sprintf("vis.tsne"))
    }
  }else{
    obj <- ssc.variableGene(obj,method=method.vgene,sd.n=sd.n,assay.name=assay.name)
    obj <- ssc.clustSubsamplingClassification(obj, assay.name=assay.name,
                                                   frac=sub.frac,
                                                   method.vgene=method.vgene,
                                                   method.reduction=method.reduction,
                                                   method.clust=method.clust,
                                                   method.classify=method.classify,
                                                   pca.npc=pca.npc,
                                                   iCor.niter=iCor.niter,
                                                   use.proj=sub.use.proj,
                                                   vis.proj=sub.vis.proj,
                                                   k.batch=k.batch)
  }
  return(obj)
}

#' Plot on tSNE map
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param gene, character; genes to be showed. (default: NULL)
#' @param columns character; columns in colData(obj) to be showd. (default: NULL)
#' @param plotDensity logical; whether plot 2D density. (default F)
#' @param colSet list; mapping iterms in the names to colors in the values. (default: list())
#' @param reduced.name character; names in the reducedDimNames. (default: "iCor.tsne")
#' @param reduced.dim integer; which dimensions of the reduced data to be used. (default: c(1,2))
#' @param out.prefix character; output prefix. (default: NULL)
#' @param p.ncol integer; number of columns in the figure layout. (default: 3)
#' @param width numeric; width of the plot, used for geneOnTSNE. (default: NA)
#' @param height numeric; height of the plot, used for geneOnTSNE. (default: NA)
#' @param base_aspect_ratio numeric; base_aspect_ratio, used for plotting metadata. (default 1.1)
#' @param peaks integer or character; index or names of the peaks. (default: NULL)
#' @param xlim integer or NULL; only draw points lie in the ragne specified by xlim and ylim (default NULL)
#' @param ylim integer or NULL; only draw points lie in the ragne specified by xlim and ylim (default NULL)
#' @param size double; points' size. If NULL, infer from number of points (default NULL)
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual theme_bw aes_string guides guide_legend
#' @importFrom cowplot save_plot plot_grid
#' @importFrom utils read.table
#' @importFrom RColorBrewer brewer.pal
#' @details If `gene` is not NULL, expression of the specified genes will be plot on the tSNE map; if columns in not
#' NULL, colData of obj with names in `columns` will be plot on the tSNE map. The tSNE map used is specified by option
#' `reduced.name` and `reduced.dim`. Both `gene` and `columns` can be non-NULL. For list `colSet`, each element define
#' a color mapping for the responding iterm in the `column`; if not specifed, automatically generated color mapping will
#' be used.
#' @export
ssc.plot.tsne <- function(obj, assay.name="exprs", gene=NULL, columns=NULL, plotDensity=F, colSet=list(),
                          reduced.name="iCor.tsne",reduced.dim=c(1,2),xlim=NULL,ylim=NULL,size=NULL,
                          out.prefix=NULL,p.ncol=3,width=NA,height=NA,base_aspect_ratio=1.1,peaks=NULL)
{
  #requireNamespace("ggplot2")
  #requireNamespace("cowplot")
  if(length(reduced.dim)!=2){ warning(sprintf("Wrong parameter, reduced.dim!!")); return(); }
  if(is.null(reducedDim(obj,reduced.name))){
    warning(sprintf("No reducedDim: %s\n",reduced.name))
    return()
  }
  dat.map <- reducedDim(obj,reduced.name)[,reduced.dim]
  if(!is.null(columns))
  {
    if(all(columns %in% colnames(colData(obj))))
    {
      if(is.list(colSet)){
        multi.p <- lapply(columns,function(cc){
          if(is.null(colSet[[cc]])){
            cc.values <- sort(unique(colData(obj)[,cc]))
            colSet[[cc]] <- structure(auto.colSet(length(cc.values),name = "Paired"),
                                      names=cc.values)
          }
          dat.plot <- data.frame(sample=rownames(dat.map),stringsAsFactors = F)
          dat.plot <- as.data.frame(cbind(dat.plot,dat.map,colData(obj)[,cc,drop=F]))
          colnames(dat.plot) <- c("sample","Dim1","Dim2",cc)
          dat.plot <- dat.plot[order(dat.plot[,cc]),]
          npts <- nrow(dat.plot)
          p <- ggplot2::ggplot(dat.plot,aes(Dim1,Dim2)) +
            geom_point(aes_string(colour=cc),size=if(is.null(size)) auto.point.size(npts)*1.1 else size)
          if(is.numeric(dat.plot[,cc])){
            p <- p + scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, "YlOrRd"))
          }else{
            p <- p + scale_colour_manual(values = colSet[[cc]])
          }
          p <- p + theme_bw() + coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE) +
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
        warning(sprintf("invalidate parameter: colSet. Please check that!"))
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
    if(is.null(names(gene))){
      names(gene) <- gene
    }
    p <- ggGeneOnTSNE(assay(obj,assay.name),
                     dat.map,
                     gene,out.prefix,p.ncol=p.ncol,xlim=xlim,ylim=ylim,size=size,width=width,height=height)
    if(is.null(out.prefix)){ print(p) }
  }
  if(plotDensity){
    if(is.null(out.prefix)){
      plot.density2D(dat.map,peaks = peaks)
    }else{
      pdf(sprintf("%s.density.pdf",out.prefix),width = 5,height = 5)
      plot.density2D(dat.map,peaks = peaks)
      dev.off()
    }
  }
}


#' Plot violin
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param gene character; genes to be showed. (default: NULL)
#' @param columns character; columns in colData(obj) to be showd. (default: NULL)
#' @param group.var character; column in the colData(obj) used for grouping. (default: "majorCluster")
#' @param out.prefix character; output prefix. (default: NULL)
#' @param p.ncol integer; number of columns in the figure layout. (default: 3)
#' @param base_aspect_ratio numeric; base_aspect_ratio, used for plotting metadata. (default 1.1)
#' @param ... parameter passed to cowplot::save_plot
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_violin scale_fill_gradient2 theme_bw theme aes_string facet_grid element_text
#' @importFrom cowplot save_plot plot_grid
#' @importFrom data.table melt
#' @details If `gene` is not NULL, violin of the genes' expression will be plot; if columns in not
#' NULL, colData of obj with names in `columns` will be plot in violin.
#' @export
ssc.plot.violin <- function(obj, assay.name="exprs", gene=NULL, columns=NULL,
                            group.var="majorCluster",
                            out.prefix=NULL,p.ncol=1,base_aspect_ratio=1.1,...)
{
  requireNamespace("ggplot2")
  requireNamespace("data.table")
  gene <- ssc.displayName2id(obj,display.name = gene)
  dat.plot <- t(assay(obj,assay.name)[gene,])
  colnames(dat.plot) <- ssc.id2displayName(obj,colnames(dat.plot))
  dat.plot.df <- data.table::data.table(sample=rownames(dat.plot),stringsAsFactors = F)
  dat.plot.df[,group.var] <- colData(obj)[,group.var]
  dat.plot.df <- cbind(dat.plot.df,dat.plot)
  dat.plot.df <- data.table::melt(dat.plot.df,id.vars=c("sample",group.var),
                                  variable.name="gene",value.name=assay.name)
  dat.plot.df.grpMean <- dat.plot.df[,lapply(.SD,mean),by=c("gene",group.var),.SDcols=assay.name]
  colnames(dat.plot.df.grpMean) <- c("gene",group.var,"meanExp")
  dat.plot.df <- dat.plot.df.grpMean[dat.plot.df,,on=c("gene",group.var)]
  dat.plot.df[meanExp<0,meanExp:=0,]
  dat.plot.df[meanExp>15,meanExp:=15,]
  head(dat.plot.df)
  p <- ggplot(dat.plot.df, aes_string(group.var, assay.name)) +
    geom_violin(scale = "width",aes(fill=meanExp),color=NA,show.legend = T) +
    scale_fill_gradient2(low = "yellow",mid = "red",high = "black",midpoint = 7.5,
                        limits=c(0,15)) +
    theme_bw(base_size = 12) +
    facet_grid(gene ~ .,switch = "y",scales = "free_y") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),strip.placement = "inside")
  if(!is.null(out.prefix)){
    cowplot::save_plot(sprintf("%s.violin.gene.pdf",out.prefix),p,
                       ncol = p.ncol,
                       base_aspect_ratio=base_aspect_ratio,...)
  }else{
    print(p)
  }
}


#' Plot pca result, such as scree plot.
#' @param obj object of \code{singleCellExperiment} class
#' @param out.prefix character; output prefix. (default: NULL)
#' @param p.ncol integer; number of columns in the figure layout. (default: 2)
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual theme_bw ylab
#' @export
ssc.plot.pca <- function(obj, out.prefix=NULL,p.ncol=2)
{
  requireNamespace("ggplot2")
  eigenv <- metadata(obj)$ssc$pca.res$eigenv.prop * 100
  dat.plot.eigenv <- data.frame(PC=seq_along(eigenv),
                                eigenv=eigenv,
                                isKneePts=as.character(seq_along(eigenv)==metadata(obj)$ssc$pca.res$kneePts),
                                stringsAsFactors = F)
  p <- ggplot2::ggplot(head(dat.plot.eigenv,n=30),mapping = aes(PC,eigenv)) +
    geom_point(aes(colour=isKneePts),show.legend=F) + ylab("Variation explained (%)") +
    scale_colour_manual(values = c("TRUE"="#E41A1C","FALSE"="#377EB8")) +
    theme_bw()
  print(p)
}

