#' @importFrom sscVis ssc.plot.silhouette
#' @export
sscVis::ssc.plot.silhouette

#' @importFrom sscVis ssc.plot.tsne
#' @export
sscVis::ssc.plot.tsne

#' @importFrom sscVis ssc.plot.violin
#' @export
sscVis::ssc.plot.violin

#' @importFrom sscVis ssc.plot.pca
#' @export
sscVis::ssc.plot.pca

#' @importFrom sscVis ssc.plot.cor
#' @export
sscVis::ssc.plot.cor

#' @importFrom sscVis ssc.plot.heatmap
#' @export
sscVis::ssc.plot.heatmap

#' @importFrom sscVis ssc.plot.GeneDensity
#' @export
sscVis::ssc.plot.GeneDensity

#' @importFrom sscVis ssc.assay.hclust
#' @export
sscVis::ssc.assay.hclust

#' @importFrom sscVis ssc.assay.zscore
#' @export
sscVis::ssc.assay.zscore

#' @importFrom sscVis ssc.order
#' @export
sscVis::ssc.order

#' @importFrom sscVis ssc.average.cell
#' @export
sscVis::ssc.average.cell

#' @importFrom sscVis ssc.displayName2id
#' @export
sscVis::ssc.displayName2id

#' @importFrom sscVis ssc.id2displayName
#' @export
sscVis::ssc.id2displayName

#' @importFrom sscVis ssc.downsample
#' @export
sscVis::ssc.downsample

#' @importFrom sscVis ssc.toLongTable
#' @export
sscVis::ssc.toLongTable

#' @importFrom sscVis ssc.moduleScore
#' @export
sscVis::ssc.moduleScore

#' @importFrom sscVis ssc.scale
#' @export
sscVis::ssc.scale

#' @importFrom sscVis ssc.build
#' @export
sscVis::ssc.build

#' @importFrom sscVis effectsize
#' @export
sscVis::effectsize

#' @importFrom sscVis directEScombi
#' @export
sscVis::directEScombi

#' @importFrom sscVis directEScombiFromLongTable
#' @export
sscVis::directEScombiFromLongTable

#' @importFrom sscVis collapseEffectSizeLong
#' @export
sscVis::collapseEffectSizeLong

#### ===================================

#' Identify variable genes
#'
#' Identify variable genes which will be used in downstream analysis. Multiple methods are available.
#'
#' @importFrom stats sd
#' @importFrom scran trendVar decomposeVar
#' @importFrom SummarizedExperiment colData rowData `colData<-` `rowData<-`
#' @importFrom S4Vectors metadata `metadata<-`
#' @importFrom utils head
#' @importFrom graphics plot points curve
#' @param obj object of SingleCellExperiment
#' @param method method to be used, can be one of "HVG.sd", "HVG.mean.sd", "HVG.trendVar". (default: "HVG.sd")
#' @param sd.n top number of genes (default 1500)
#' @param mean.thre numeric; threshold for mean, used in trendVar method (default 0.1)
#' @param fdr.thre numeric; threshold for fdr, used in trendVar method (default 0.001)
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
ssc.variableGene <- function(obj,method="HVG.sd",sd.n=1500,mean.thre=0.1,fdr.thre=0.001,assay.name="exprs",
                             var.block=NULL,
                             reuse=F,out.prefix=NULL)
{
  if(method=="HVG.sd")
  {
    ###if(!reuse || is.null(metadata(obj)$ssc[["variable.gene"]][["sd"]])){
    if(!reuse || !("HVG.sd" %in% colnames(rowData(obj))) ){
      rowData(obj)[,"row.sd"] <- apply(assay(obj,assay.name),1,sd)
      rowData(obj)[,"HVG.sd"]  <- rownames(obj) %in% names(head(sort(rowData(obj)[,"row.sd"],decreasing = T),n=sd.n))
    }
  }else if(method=="HVG.mean.sd")
  {
    if(!reuse || !("HVG.mean.sd" %in% colnames(rowData(obj))) ){
      rowData(obj)[,"row.sd"] <- apply(assay(obj,assay.name),1,sd)
      rowData(obj)[,"row.mean"] <- apply(assay(obj,assay.name),1,mean)
      f.gene <- rowData(obj)[,"row.sd"]>1 & rowData(obj)[,"row.mean"]>1
      info.gene.sd <- rowData(obj)[f.gene,"row.sd"]
      info.gene.sd <- sort(info.gene.sd,decreasing = T)
      rowData(obj)[,"HVG.mean.sd"] <- rownames(obj) %in% names(head(info.gene.sd,n=sd.n))
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
    f.var <- var.out$FDR<fdr.thre & var.out$mean>mean.thre
    rowData(obj)[["trendVar"]] <- rownames(obj) %in% rownames(var.out)[f.var]
    metadata(obj)$ssc[["variable.gene"]][["trendVar.detail"]] <- var.out
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


#######################################


#' Reduce dimension by various methods
#' @importFrom SingleCellExperiment reducedDim `reducedDim<-` reducedDimNames
#' @importFrom S4Vectors metadata `metadata<-`
#' @importFrom stats cor prcomp
#' @importFrom BiocGenerics t
#' @importFrom uwot umap
#' @importFrom utils head
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param method character; method to be used for dimension reduction, should be one of (pca, tsne, iCor). (default: "iCor")
#' @param method.vgene character; method to identify variable genes. (default: sd)
#' @param method.tsne character; method to run tsne, one of "Rtsne", "FIt-SNE". (default: "Rtsne")
#' @param pca.npc integer; number of pc be used. Only for reduction method "pca". (default: NULL)
#' @param tSNE.usePCA logical; whether use PCA before tSNE. Only for reduction method "tsne". (default: T)
#' @param tSNE.perplexity logical; perplexity parameter. Used in all Rtsne() calling. (default: 30)
#' @param autoTSNE logical; Wheter generate automatically a tSNE map when reduction method is "pca" or "iCor". (default: T)
#' @param dim.name character; store the reduced data under the name in the obj's reducedDim SimpleList. If i
#' it is NULL, infer from method. (default: NULL)
#' @param iCor.niter integer; number of iteration of calculating the correlation. Used in reduction method "iCor". (default: 1)
#' @param iCor.method character; method to calculate correlation between samples,
#' should be one of "spearman" and "pearson". (default "spearman")
#' @param reuse logical; don't calculate if the query is already available. (default: F)
#' @param seed integer; seed of random number generation. (default: NULL)
#' @param out.prefix character; output prefix (default: NULL)
#' @param ... parameters passed to RD methods
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
                          method.vgene="HVG.sd",
                          method.tsne="Rtsne",
                          pca.npc=NULL,
                          tSNE.usePCA=T,
                          tSNE.perplexity=30,
                          autoTSNE=T,
                          dim.name=NULL,
                          iCor.niter=1,iCor.method="spearman",
                          reuse=F,ncore=NULL,seed=NULL,out.prefix=NULL,...)
{
  row.sd <- apply(assay(obj,assay.name),1,sd)
  col.sd <- apply(assay(obj,assay.name),2,sd)
  col.zero <- which(col.sd==0)
  if(any(is.na(row.sd))){ warning(sprintf("expression data contians na value\n")) }
  if(length(col.zero>0)){
    warning(sprintf("expression data contains colum(s) with sd equal to 0:\n%s\n",
                    paste(head(colnames(obj)[col.zero]),collapse = ",")))
  }
  obj <- obj[!is.na(row.sd) & row.sd>0,]
  if(!method.vgene %in% colnames(rowData(obj)) ){
    stop(sprintf("No variable genes identified by method %s !",method.vgene))
  }
  if(is.null(dim.name)){ dim.name <- method }
  vgene <- rowData(obj)[[method.vgene]]
  if(!reuse || !(dim.name %in% reducedDimNames(obj)) ){
    if(is.null(seed)){
      myseed <- as.integer(Sys.time())
    }else{
      myseed <- seed
    }
    loginfo(sprintf("set.seed(%s) for ssc.reduceDim\n",as.character(myseed)))
    set.seed(myseed)

    if(method=="pca"){
        pca.res <- prcomp(BiocGenerics::t(assay(obj[vgene,],assay.name)))
        ### find elbow point and get number of components to be used
        pca.res$eigenv.prop <- (pca.res$sdev^2)/sum(pca.res$sdev^2)
        pca.res$eigengap <- pca.res$eigenv.prop[-length(pca.res$eigenv.prop)]-pca.res$eigenv.prop[-1]
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
        if(autoTSNE){
            reducedDim(obj,sprintf("%s.tsne",dim.name)) <-
                run.tSNE(proj_data,tSNE.usePCA=F,tSNE.perplexity,
                         out.prefix=out.prefix,n.cores=ncore,method=method.tsne,...)
        }
    }else if(method=="tsne"){
      proj_data <- run.tSNE(BiocGenerics::t(assay(obj[vgene,],assay.name)),tSNE.usePCA,tSNE.perplexity,
                            out.prefix=out.prefix,n.cores=ncore,method=method.tsne,...)
    }else if(method=="umap"){
        proj_data <- umap(BiocGenerics::t(assay(obj[vgene,],assay.name)), init = "spca",pca=50,n_threads = ncore)
    }else if(method=="iCor"){
      proj_data <- as.matrix(assay(obj[vgene,],assay.name))
      while(iCor.niter>0){
        ##proj_data <- cor(proj_data,method=iCor.method)
        proj_data <- cor.BLAS(t(proj_data),method=iCor.method,nthreads = ncore)
        iCor.niter <- iCor.niter-1
      }
      if(autoTSNE) { reducedDim(obj,sprintf("%s.tsne",dim.name)) <-
          run.tSNE(proj_data,tSNE.usePCA=F,tSNE.perplexity,out.prefix=out.prefix,
                   n.cores=ncore,method=method.tsne,...)
      }
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
#' @importFrom leiden leiden
#' @importFrom SingleCellExperiment reducedDimNames `reducedDim<-` reducedDim
#' @importFrom SummarizedExperiment rowData `rowData<-` colData `colData<-`
#' @importFrom S4Vectors metadata `metadata<-`
#' @importFrom stats quantile
#' @importFrom graphics plot points
#' @importFrom grDevices pdf png dev.off
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param method.reduction character; which dimention reduction method to be used, should be one of
#' "iCor", "pca" and "none". (default: "iCor")
#' @param method character; clustering method to be used, should be one of "kmeans", "hclust", "dynamicTreeCut", "SNN", "dpclust", "adpclust" and "SC3". (default: "kmeans")
#' @param k.batch integer; number of clusters to be evaluated. (default: 2:6)
#' @param method.vgene character; variable gene identification method used. (default: "HVG.sd")
#' @param SNN.k integer; number of shared NN. (default: 10)
#' @param SNN.method character; cluster method applied on SNN, one of "greedy", "eigen", "infomap",
#' "prop", "louvain", "optimal", "spinglass", "walktrap", "betweenness", "leiden". (default: "eigen")
#' @param SC3.biology logical, SC3 parameter, whether calcualte biology. (default: T)
#' @param SC3.markerplot.width integer, SC3 parameter, with of the marker plot (default: 15)
#' @param dpclust.rho numberic; cuttoff of rho, if it is NULL, infer frome the data (default: NULL)
#' @param dpclust.delta numberic; cuttoff of delta, if it is NULL, infer frome the data (default: NULL)
#' @param out.prefix character; output prefix, if not NULL, some plots of intermediate result will be produced. (default: NULL)
#' @param parlist list; if not NULL, use th parameters in it. (default: NULL)
#' @param seed integer; seed of random number generation. (default: NULL)
#' @param ncore integer; nuber of CPU cores to use. if NULL, automatically detect the number. (default: NULL)
#' @param ... parameters pass to clustering methods
#' @details If no dimension reduction performed or method is "none", expression data of variable genes,
#' which can be speficed by method.vgene, will be used for clustering. Otherwise, the reduced data specified by
#' method.reduction will be used. The cluster label will stored in the colData of the object
#' of \code{singleCellExperiment} class, with colname in the format of \{method.reduction\}.\{method\}k\{k\}
#' where \{k\} get value(s) from k.batch.
#' @return an object of \code{SingleCellExperiment} class with cluster labels added.
#' @export
ssc.clust <- function(obj, assay.name="exprs", method.reduction="iCor",
                      method="kmeans", k.batch=2:6,
                      method.vgene="HVG.sd",
                      SNN.k=10,SNN.method="eigen",
                      SC3.biology=T,SC3.markerplot.width=15,
                      dpclust.rho=NULL,dpclust.delta=NULL,
                      parlist=NULL,
                      out.prefix=NULL,seed=NULL,ncore=NULL,...)
{
  clust.res <- NULL
  res.list <- list()
  ### check transformed data
  if(method.reduction=="none"){
    warning(sprintf("The dimention reduction should be performed!"))
    vgene <- rowData(obj)[[method.vgene]]
    dat.transformed <- t(assay(obj[vgene,],assay.name))
  }else if( (!method.reduction %in% reducedDimNames(obj)) ){
    dat.transformed <- NULL
  }else{
    dat.transformed <- reducedDim(obj,method.reduction)
  }
  if(is.null(dat.transformed) && method!="SC3"){
    warning("dat.transformed is null !!")
    metadata(obj)$ssc$clust.res[[method]] <- NULL
    if(method %in% c("adpclust","dpclust","SNN","dynamicTreeCut")){
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
    colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%02d",clust.res$clusters)
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
    colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%02d",clust.res$clusters)
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
    }else if(SNN.method=="leiden"){
      clust.res <- leiden::leiden((snn.gr),seed=myseed,...)
	  clust.res <- list(membership=clust.res)
	}
    metadata(obj)$ssc$clust.res[["snn.gr"]] <- snn.gr
    loginfo(sprintf("cluster using method %s finished\n",SNN.method))
    k <- "auto"
	colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%02d",clust.res$membership)
    res.list[[as.character(k)]] <- clust.res
  }else if(method=="SC3"){
    ##vgene <- metadata(obj)$ssc[["variable.gene"]][[method.vgene]]
    vgene <- rowData(obj)[[method.vgene]]
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
  }else if(method=="dynamicTreeCut"){
    k <- "auto"
    obj.distM <- stats::dist(dat.transformed)
    obj.hclust <- stats::hclust(obj.distM,method="ward.D2")
    cluster.label <- dynamicTreeCut::cutreeDynamic(obj.hclust,distM=as.matrix(obj.distM),
                                                   method = "hybrid", ...)
    clust.res <- list("hclust"=obj.hclust,"dist"=obj.distM,"cluster"=cluster.label)
    colData(obj)[,sprintf("%s.%s.k%s",method.reduction,method,k)] <- sprintf("C%02d",clust.res$cluster)
    res.list[[as.character(k)]] <- clust.res
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
            clust.res <- factoextra::eclust(dat.transformed, "hclust",
                                            k = if(k==0) NULL else k, graph = FALSE,...)
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
#' @importFrom SingleCellExperiment reducedDims `reducedDims<-` `reducedDim<-`
#' @importFrom SummarizedExperiment rowData colData `rowData<-` `colData<-`
#' @importFrom S4Vectors SimpleList metadata `metadata<-`
#' @importFrom stats cor
#' @importFrom plyr llply
#' @importFrom MASS ginv
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param frac numeric; subsample to frac of original samples. (default: 0.4)
#' @param method.vgene character; variable gene identification method used. (default: "HVG.sd")
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
                                               method.vgene="HVG.sd",
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
  #metadata(obj.train)$ssc[["variable.gene"]][[method.vgene]] <- intersect(metadata(obj.train)$ssc[["variable.gene"]][[method.vgene]],
  #                                                                      rownames(obj.train))
  #metadata(obj.pred)$ssc[["variable.gene"]][[method.vgene]] <- intersect(metadata(obj.pred)$ssc[["variable.gene"]][[method.vgene]],
  #                                                                       rownames(obj.pred))

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
    #vgene <- metadata(obj.pred)$ssc$variable.gene[[method.vgene]]
    vgene <- rowData(obj.pred)[[method.vgene]]

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
#' @param fdr.thre numeric; threshold for fdr, used in trendVar method (default 0.001)
#' @param var.block character; specify the uninteresting factors by formula. E.g. "~patient".
#' used in trendVar method (default NULL)
#' @param sd.n integer; top number of genes as variable genes (default 1500)
#' @param de.n integer; number of differential genes used for refined geneset for another run of clustering (default 1500)
#' @param method.reduction character; which dimention reduction method to be used, should be one of
#' "iCor", "pca", and "none". (default: "iCor")
#' @param method.clust character; clustering method to be used, should be one of "kmeans", "hclust", "SNN", "adpclust" and "SC3. (default: "kmeans")
#' @param method.classify character; method used for classification, one of "knn" and "RF". (default: "knn")
#' @param method.tsne character; method to run tsne, one of "Rtsne", "FIt-SNE". (default: "Rtsne")
#' @param pca.npc integer; number of pc be used. Only for reduction method "pca". (default: NULL)
#' @param iCor.niter integer; number of iteration of calculating the correlation. Used in reduction method "iCor". (default: 1)
#' @param iCor.method character; correlation method, one of "spearman", "pearson" (default: "spearman")
#' @param tSNE.perplexity double, perplexity parameter of tSNE. (default: 30)
#' @param subsampling logical; whether cluster using the subsampling->cluster->classification method. (default: F)
#' @param sub.frac numeric; subsample to frac of original samples. (default: 0.4)
#' @param sub.use.proj logical; whether use the projected data for classification. (default: T)
#' @param sub.vis.proj logical; whether get low dimensional representation for visualization, only used in downsample mode. (default: F)
#' @param k.batch integer; number of clusters to be evaluated. (default: 2:6)
#' @param refineGene logical; whether perform second round demension reduction and clustering pipeline using the differential
#' genes found by the first round cluster result. (default: F)
#' @param HSD.FC.THRESHOLD numeric; threshold for log2FoldChange, used in findDEGenesByAOV (default 1)
#' @param nIter integer; number of iterative clustering in sub-cluster. (default: 1)
#' @param do.DE logical; perform DE analysis when clustering finished. (default: F)
#' @param out.prefix character; output prefix, if not NULL, some plots of intermediate result will be produced. (default: NULL)
#' @param parfile character; parameter files, if not NULL, will use the settings. must contain a list named
#' `parlist`. (default: NULL)
#' @param ncore integer; nuber of CPU cores to use. if NULL, automatically detect the number. (default: NULL)
#' @param reuse logical; don't calculate if the query is already available. (default: F)
#' @param seed integer; seed of random number generation. (default: NULL)
#' @param ... parameters pass to clustering methods
#' @importFrom SummarizedExperiment colData rowData `colData<-` `rowData<-`
#' @importFrom S4Vectors metadata `metadata<-`
#' @importFrom utils head
#' @importFrom stats as.formula
#' @details run the pipeline of variable gene identification, dimension reduction, clustering.
#' @seealso \code{\link{ssc.variableGene}} for variable genes' identification, \code{\link{ssc.reduceDim}}
#' for dimension reduction, \code{\link{ssc.clust}} for clustering using all data
#' and \code{\link{ssc.clustSubsamplingClassification}} for clustering with subsampling.
#' @return an object of \code{SingleCellExperiment} class with cluster labels added.
#' @export
ssc.run <- function(obj, assay.name="exprs",
                    method.vgene="HVG.sd",
                    sd.n=1500,
                    mean.thre=0.1,
                    fdr.thre=0.001,
                    var.block=NULL,
                    method.reduction="iCor",
                    method.clust="kmeans",
                    method.classify="knn",
                    method.tsne="Rtsne",
                    pca.npc=NULL,
                    iCor.niter=1,
                    iCor.method="spearman",
                    tSNE.perplexity=30,
                    subsampling=F,
                    sub.frac=0.4,
                    sub.use.proj=T,
                    sub.vis.proj=F,
                    k.batch=2:6,
                    refineGene=F,
                    de.n=1500,
                    HSD.FC.THRESHOLD=1,
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
      obj <- ssc.variableGene(obj,method=method.vgene,sd.n=sd.n,mean.thre = mean.thre,fdr.thre=fdr.thre,
                              assay.name=assay.name,reuse = reuse,var.block = var.block,
                              out.prefix = sprintf("%s.%s",out.prefix,rid))
      loginfo(sprintf("reduce dimensions ... (%s)",rid))
      obj <- ssc.reduceDim(obj,assay.name=assay.name,
                           method=method.reduction,
                           pca.npc = pca.npc,
                           iCor.niter = iCor.niter,
                           iCor.method = iCor.method,
                           method.vgene=method.vgene,
                           method.tsne=method.tsne,tSNE.perplexity=tSNE.perplexity,
                           ncore = ncore,
                           seed = seed,
                           reuse = reuse)
      loginfo(sprintf("clustering ... (%s)",rid))
      obj <- ssc.clust(obj, assay.name=assay.name,
                       method.reduction=if(method.clust %in% c("adpclust","dpclust") && method.reduction!="umap") sprintf("%s.tsne",method.reduction) else method.reduction,
                       method=method.clust, k.batch=k.batch,
                       out.prefix = if(is.null(out.prefix)) NULL else sprintf("%s.%s",out.prefix,rid),
                       seed = seed,
                       method.vgene=method.vgene, parlist = parlist.rid, ...)

      .xlabel <- NULL
      if(method.clust %in% c("adpclust","dpclust") && method.reduction!="umap"){
        .xlabel <- sprintf("%s.tsne.%s.kauto",method.reduction,method.clust)
      ##}else if(method.clust=="SNN"){
      }else{
        .xlabel <- sprintf("%s.%s.kauto",method.reduction,method.clust)
      }
      if(!is.null(out.prefix) && !is.null(.xlabel)){
        ssc.plot.tsne(obj,columns = c(.xlabel),
                      reduced.name = if(method.clust %in% c("adpclust","dpclust") && method.reduction!="umap") sprintf("%s.tsne",method.reduction) else method.reduction,
                      out.prefix = sprintf("%s.%s",out.prefix,rid),
                      base_aspect_ratio = 1.4)
        if(method.reduction=="pca"){
            pdf(sprintf("%s.%s.pca.00.pdf",out.prefix,rid),width = 8,height = 6)
            ssc.plot.pca(obj)
            dev.off()
        }
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
                                     HSD.FC.THRESHOLD = HSD.FC.THRESHOLD,
                                     gid.mapping = rowData(obj)[,"display.name"])
          if(!is.null(de.out) && nrow(de.out$aov.out.sig)>30)
          {
            metadata(obj)$ssc[["de.res"]][[rid]] <- de.out
            #metadata(obj)$ssc[["variable.gene"]][["refine.de"]] <- head(de.out$aov.out.sig$geneID,n=de.n)
            rowData(obj)[,"HVG.refine.de"]  <- rownames(obj) %in% head(de.out$aov.out.sig$geneID,n=de.n)
            loginfo(sprintf("reduce dimensions using DE genes ... (%s)",rid))
            obj <- ssc.reduceDim(obj,assay.name=assay.name,
                         method=method.reduction,
                         method.tsne=method.tsne,
                         pca.npc = pca.npc,
                         iCor.niter = iCor.niter,
                         iCor.method = iCor.method,
                         tSNE.perplexity=tSNE.perplexity,
                         ncore = ncore,
                         seed = seed,
                         dim.name = sprintf("de.%s",method.reduction),
                         method.vgene="HVG.refine.de",reuse = reuse)
            loginfo(sprintf("clust using DE genes ... (%s)",rid))
            obj <- ssc.clust(obj, assay.name=assay.name,
                             method.reduction=if(method.clust %in% c("adpclust","dpclust") && method.reduction!="umap") sprintf("de.%s.tsne",method.reduction) else sprintf("de.%s",method.reduction),
                             method=method.clust, k.batch=k.batch,
                             seed = seed,
                             out.prefix = if(is.null(out.prefix)) NULL else sprintf("%s.%s.refineG",out.prefix,rid),
                             method.vgene="HVG.refine.de", parlist = parlist.rid, ...)

            if(method.clust %in% c("adpclust","dpclust") && method.reduction!="umap"){
              colData(obj)[,.xlabel] <- colData(obj)[,sprintf("de.%s.tsne.%s.kauto",method.reduction,method.clust)]
              colData(obj)[,sprintf("de.%s.tsne.%s.kauto",method.reduction,method.clust)] <- NULL
            }else{
              colData(obj)[,.xlabel] <- colData(obj)[,sprintf("de.%s.%s.kauto",method.reduction,method.clust)]
              colData(obj)[,sprintf("de.%s.%s.kauto",method.reduction,method.clust)] <- NULL
            }
            if(!is.null(out.prefix) && !is.null(.xlabel)){
              ssc.plot.tsne(obj,columns = c(.xlabel),
                            reduced.name = if(method.clust %in% c("adpclust","dpclust") && method.reduction!="umap") sprintf("de.%s.tsne",method.reduction) else sprintf("de.%s",method.reduction),
                            out.prefix = sprintf("%s.%s.refineG",out.prefix,rid),
                            base_aspect_ratio = 1.4)
              if(method.reduction=="pca"){
                pdf(sprintf("%s.%s.refineG.pca.00.pdf",out.prefix,rid),width = 8,height = 6)
                ssc.plot.pca(obj)
                dev.off()
              }
            }
          }else{
            warning("The number of DE genes is less than 30, NO second round clustering using DE genes will be performed!!")
          }
        }
      }
      if(level<nIter){
        clustLabel <- if(method.clust %in% c("adpclust","dpclust") && method.reduction!="umap") sprintf("%s.tsne.%s.kauto",method.reduction,method.clust) else sprintf("%s.%s.kauto",method.reduction,method.clust)
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
      if(method.clust %in% c("adpclust","dpclust") && method.reduction!="umap"){
        .xlabel <- colData(obj)[,sprintf("%s.tsne.%s.kauto",method.reduction,method.clust)]
      }else{
        .xlabel <- colData(obj)[,sprintf("%s.%s.kauto",method.reduction,method.clust)]
      }
      de.out <- findDEGenesByAOV(xdata = assay(obj,assay.name),
                                 xlabel = .xlabel,
                                 HSD.FC.THRESHOLD = HSD.FC.THRESHOLD,
                                 gid.mapping = rowData(obj)[,"display.name"])
      metadata(obj)$ssc[["de.res"]][["L1C1"]] <- de.out
      #metadata(obj)$ssc[["variable.gene"]][["de"]] <- head(de.out$aov.out.sig$geneID,n=sd.n)
      rowData(obj)[,"HVG.de"]  <- rownames(obj) %in% head(de.out$aov.out.sig$geneID,n=sd.n)
      ### for general visualization
      obj <- ssc.reduceDim(obj,assay.name=assay.name, method="tsne",method.tsne=method.tsne,
                           tSNE.perplexity=tSNE.perplexity,
                           method.vgene="HVG.de",dim.name = sprintf("vis.tsne"))
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




#' identify marker genes of each cluster
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param group.var character; column in the colData(obj) used for grouping. (default: "majorCluster")
#' @param batch character; covariate. (default: NULL)
#' @param assay.bin character; binarized expression assay (default: NULL)
#' @param out.prefix character; output prefix. (default: NULL)
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
#' @param do.plot logical; whether plot. (default: TRUE)
#' @param F.FDR.THRESHOLD numeric; threshold of the adjusted p value of F-test. (default: 0.01)
#' @param pairwise.P.THRESHOLD numeric; threshold of the adjusted p value of HSD-test (default: 0.01)
#' @param pairwise.FC.THRESHOLD numeric; threshold of the absoute diff of HSD-test (default: 1)
#' @param use.Kruskal logical; whether use Kruskal test for ranking genes (default: FALSE)
#' @param method.Max character; method to find highest group, one of "mean", "median", "rank.mean" (default: mean)
#' @param do.force logical; . (default: FALSE)
#' @param verbose logical; whether output all genes' result. (default: FALSE)
#' @param ... parameters passed to findDEGenesByAOV
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply
#' @importFrom dplyr inner_join
#' @importFrom utils write.table
#' @importFrom SummarizedExperiment rowData `rowData<-` colData `colData<-`
#' @details identify marker genes based on aov and AUC.
#' @export
ssc.clusterMarkerGene <- function(obj, assay.name="exprs", ncell.downsample=NULL,
                                  group.var="majorCluster",batch=NULL,
                                  assay.bin=NULL, out.prefix=NULL,n.cores=NULL, do.plot=T,
                                  F.FDR.THRESHOLD=0.01,pairwise.P.THRESHOLD=0.01,pairwise.FC.THRESHOLD=1,
				  use.Kruskal=F,method.Max="mean",
				  do.force=F, verbose=F,...)
{
#    requireNamespace("doParallel")
    requireNamespace("plyr")
    requireNamespace("dplyr")

    #### downsample cells
    if(!is.null(ncell.downsample)){
        clust <- colData(obj)[,group.var]
        names(clust) <- colnames(obj)
        grp.list <- unique(clust)
        f.cell <- unlist(sapply(grp.list,function(x){
                             x <- names(clust[clust==x])
                             sample(x,min(length(x),ncell.downsample)) }))
        obj <- obj[,f.cell]
    }

    clust <- colData(obj)[,group.var]
    batchV <- NULL
    if(!is.null(batch)){
        batchV <- colData(obj)[,batch]
    }
    if(!do.force && (length(unique(clust))<2 || is.null(obj) || !all(table(clust) > 1))){
        cat("WARN: clusters<2 or no obj provided or not all clusters have more than 1 samples\n")
        return(NULL)
    }
    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = n.cores)

    if(is.null(names(rowData(obj)$display.name))){ names(rowData(obj)$display.name) <- row.names(obj) }
    gid.mapping <- rowData(obj)$display.name
    if(is.null(names(gid.mapping))){ names(gid.mapping) <- row.names(obj) }
    dat.to.test <- as.matrix(assay(obj,assay.name))
    #### filter genes
    exp.frac <- NULL
    minCell <- 5
    if(!is.null(assay.bin) && assay.bin %in% assayNames(obj)){
        #.f.gid <- intersect(rownames(exp.bin),rownames(dat.to.test))
        exp.frac <- expressedFraction(as.matrix(assay(obj,assay.bin)),clust,n.cores)
        #print(str(exp.frac))
        #print(head(exp.frac))
#        f.bin <- apply(exp.frac,1,function(x){ nE <- sum(x>0.05); nNE <- sum(x<0.95); return(nE > 0 && nNE > 0) })
#        dat.to.test <- dat.to.test[f.bin,]
        f.notAllZero <- apply(dat.to.test,1,
                              function(x){ nE <- sum(x>0); return( nE > minCell & nE/length(x) > 0.01 ) })
        dat.to.test <- dat.to.test[f.notAllZero,]
    }else{
        ### currently the same with using assay.bin
        f.notAllZero <- apply(dat.to.test,1,
                              function(x){ nE <- sum(x>0); return( nE > minCell & nE/length(x) > 0.01 ) })
        dat.to.test <- dat.to.test[f.notAllZero,]
    }
    ### AUC
    .gene.table <- plyr::ldply(rownames(dat.to.test),function(x){
	sscClust:::getAUC(dat.to.test[x,],clust,use.rank = (method.Max=="rank.mean"), method.Max=method.Max,geneID=x)
    },.progress = "none",.parallel=T)
    #if(is.character(.gene.table$AUC)){ .gene.table$AUC <- as.numeric(.gene.table$AUC) }
    #if(is.character(.gene.table$score.p.value)){ .gene.table$score.p.value <- as.numeric(.gene.table$score.p.value) }
    #if(is.numeric(.gene.table$cluster)){ .gene.table$cluster <- sprintf("C%s",.gene.table$cluster) }
    .gene.table$score.q.value <- 1
    .gene.table <- merge(data.frame(geneID=rownames(dat.to.test),
                                    geneSymbol=gid.mapping[rownames(dat.to.test)],
                                    stringsAsFactors = F),
                         .gene.table)
    #### diff genes
    .aov.res <- findDEGenesByAOV(xdata = dat.to.test,
                                 xlabel =if(is.numeric(clust)) sprintf("C%s",clust) else clust,
                                 batch=batchV,
                                 out.prefix=out.prefix,
                                 F.FDR.THRESHOLD=F.FDR.THRESHOLD,
                                 HSD.FDR.THRESHOLD=pairwise.P.THRESHOLD,
                                 HSD.FC.THRESHOLD=pairwise.FC.THRESHOLD,
				 use.Kruskal=use.Kruskal,
                                 verbose=verbose,n.cores=n.cores, gid.mapping=gid.mapping,...)
    if(verbose){
        .gene.table <- dplyr::inner_join(x = .gene.table,y = .aov.res$aov.out)
    }else{
        .gene.table <- dplyr::inner_join(x = .gene.table,y = .aov.res$aov.out.sig)
    }
    ### average expression of each cluster
    avg.exp <- plyr::ldply(as.character(.gene.table$geneID),function(v){
                .res <- aggregate(dat.to.test[v,],by=list(clust),FUN=mean)
                .res.2 <- aggregate(dat.to.test[v,],by=list(clust),FUN=sd)
                structure(c(.res[,2],.res.2[,2]),names=c(sprintf("avg.%s",.res[,1]),sprintf("sd.%s",.res[,1])))
			},.progress = "none",.parallel=T)
    ##rownames(avg.exp) <- rownames(dat.to.test[as.character(.gene.table$geneID),,drop=F])
    rownames(avg.exp) <- as.character(.gene.table$geneID)
    ####colnames(avg.exp) <- sprintf("avg.%s",colnames(avg.exp))
    avg.exp.df <- data.frame(geneID=rownames(avg.exp),
                             geneSymbol=gid.mapping[rownames(avg.exp)],
                             stringsAsFactors = F)
    avg.exp.df <- cbind(avg.exp.df,avg.exp)
    .gene.table <- dplyr::inner_join(x=.gene.table,y=avg.exp.df)
    ### binarized expression fraction
    if(!is.null(exp.frac)){
        exp.frac.df <- data.frame(geneID=rownames(exp.frac),
                                  geneSymbol=gid.mapping[rownames(exp.frac)],
                                  stringsAsFactors = F)
        exp.frac.df <- cbind(exp.frac.df,exp.frac)
        .gene.table <- dplyr::inner_join(x = .gene.table,y = exp.frac.df)
    }
    rownames(.gene.table) <- .gene.table$geneID
    if(verbose && !is.null(out.prefix)){
        .gene.table <- .gene.table[order(.gene.table$F,decreasing = T),]
        write.table(.gene.table, file = sprintf("%s.geneTable.all.txt",out.prefix),
                    row.names = F,quote = F,sep = "\t")
    }
    ### final .gene.table always contain only thoase show significant differential expression
    f.isSig <- intersect(.aov.res$aov.out.sig$geneID,rownames(.gene.table))
    .gene.table <- .gene.table[f.isSig,]
    .gene.table$score.q.value <- p.adjust(.gene.table$score.p.value,method = "BH")
    order.gene <- order(.gene.table$cluster,-.gene.table$AUC)
    .gene.table <- .gene.table[order.gene,,drop=F]

    if(!is.null(out.prefix)){
        #rownames(.gene.table) <- .gene.table$geneID
        ### save .txt file
        write.table(.gene.table, file = sprintf("%s.markerGene.all.txt",out.prefix),
                    row.names = F,quote = F,sep = "\t")
        #write.table(subset(.gene.table,score.q.value<0.01),
        #            file = sprintf("%s.markerGene.q01.txt",out.prefix),
        #            row.names = F,quote = F,sep = "\t")
    }
    if(do.plot)
    {

    }
    return(list(gene.table=.gene.table,aov.res=.aov.res))
}


#' identify differential genes of each cluster (comparing the cluster with all others), using limma
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param group.var character; column in the colData(obj) used for grouping. (default: "majorCluster")
#' @param group.list character; DEG of groups to calculate. If NULL, all groups. (default: "NULL")
#' @param group.mode character; One of "multi", "multiAsTwo" (default: "multi")
#' @param batch character; covariate. (default: NULL)
#' @param out.prefix character; output prefix. (default: NULL)
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
#' @param do.plot logical; whether plot. (default: TRUE)
#' @param T.fdr numeric; threshold of the adjusted p value of moderated t-test (default: 0.05)
#' @param T.logFC numeric; threshold of the absoute diff (default: 1)
#' @param T.expr numeric; threshold for binarizing exprs (default: 0.3)
#' @param T.bin.useZ logical; wheter use the z-score version of assay.namme for binarizing exprs (default: T)
#' @param verbose integer; verbose level. (default: 0)
#' @param do.force logical; . (default: FALSE)
#' @param method character; . (default: "limma")
#' @importFrom SummarizedExperiment rowData colData `rowData<-` `colData<-`
#' @importFrom utils write.table
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom limma lmFit eBayes topTable
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply llply
#' @details identify differential genes using limma
#' @export
ssc.DEGene.limma <- function(obj, assay.name="exprs", ncell.downsample=NULL,
                                  group.var="majorCluster",group.list=NULL,group.mode="multi",batch=NULL,
                                  out.prefix=NULL,n.cores=NULL, do.plot=T,
                                  T.fdr=0.01,T.logFC=1,T.expr=0.3,T.bin.useZ=T,
                                  verbose=0,do.force=F,method="limma")
{
#    requireNamespace("doParallel")
    requireNamespace("plyr")
    requireNamespace("dplyr")

    if(!is.null(ncell.downsample) && group.mode=="multi"){
        obj <- ssc.downsample(obj, ncell.downsample=ncell.downsample, group.var=group.var,rn.seed=9999)
    }

    clust <- colData(obj)[,group.var]
    if(is.null(group.list)){ group.list <- unique(clust) }
    batchV <- NULL
    if(!is.null(batch)  && length(unique(colData(obj)[,batch])) > 1 ){
        batchV <- colData(obj)[,batch]
    }
    stat.clust <- table(clust)
    if(length(unique(clust))<2 || is.null(obj) ){
        cat("WARN: clusters<2 or no obj provided or not all clusters have more than 1 samples\n")
        return(NULL)
    }
    if(!do.force && !all(stat.clust[group.list]>=2)){
        cat("WARN: not all clusters specified have more than 2 samples\n")
        return(NULL)
    }

    if(is.null(names(rowData(obj)$display.name))){ names(rowData(obj)$display.name) <- row.names(obj) }
    gid.mapping <- rowData(obj)$display.name
    if(is.null(names(gid.mapping))){ names(gid.mapping) <- row.names(obj) }

    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = n.cores)
    out <- llply(group.list,function(x){
        xlabel <- clust
	xgroup <- x
	if(group.mode=="multiAsTwo" || method!="limma"){
	    xlabel <- as.character(xlabel)
	    xlabel[xlabel!=x] <- "_control"
	    xlabel[xlabel==x] <- "_case"
	    xlabel <- factor(xlabel,levels=c("_control","_case"))
	    xgroup <- sprintf("_case:%s",x)
	}
	if(method=="limma")
	{
	    obj.out <- run.limma.matrix(assay(obj,assay.name),xlabel,batch=batchV,
					  out.prefix=if(is.null(out.prefix)) NULL else sprintf("%s.%s",out.prefix,x),
					  group=xgroup,
					  rn.seed=9999,
					  ncell.downsample=if(group.mode=="multi") NULL else ncell.downsample,
					  T.fdr=T.fdr,T.logFC=T.logFC,T.expr=T.expr,T.bin.useZ=T.bin.useZ,
					  verbose=verbose,n.cores=1,
					  gid.mapping=gid.mapping, do.voom=F)
	}else{
	    obj.out <- run.DE.matrix(assay(obj,assay.name),xlabel,batch=batchV,
					  out.prefix=if(is.null(out.prefix)) NULL else sprintf("%s.%s",out.prefix,x),
					  group=xgroup,
					  rn.seed=9999,
					  ncell.downsample=if(group.mode=="multi") NULL else ncell.downsample,
					  T.fdr=T.fdr,T.logFC=T.logFC,T.expr=T.expr,T.bin.useZ=T.bin.useZ,
					  verbose=verbose,n.cores=1,
					  gid.mapping=gid.mapping, do.voom=F,method=method)
	}
	return(obj.out)

    },.parallel=T)
    names(out) <- group.list

    all.table <- data.table(ldply(group.list,function(x){ out[[as.character(x)]]$all }))
    sig.table <- data.table(ldply(group.list,function(x){ out[[as.character(x)]]$sig }))

    fit.list <- NULL

    if(verbose>1 & method=="limma"){
	fit.list <- llply(group.list,function(x){ out[[as.character(x)]]$fit })
	names(fit.list) <- group.list
    }

    if(verbose>2)
    {
	res.aov <- findDEGenesByAOV(as.matrix(assay(obj,assay.name)),clust,batch=batchV, out.prefix=NULL,
								n.cores=n.cores, gid.mapping=gid.mapping)
	res.aov$aov.out$F.rank <- rank(-res.aov$aov.out$F)/nrow(res.aov$aov.out)
	all.table <- merge(all.table,res.aov$aov.out[,c("geneID","F","F.pvalue","F.adjp","F.rank")],by="geneID")
	sig.table <- merge(sig.table,res.aov$aov.out[,c("geneID","F","F.pvalue","F.adjp","F.rank")],by="geneID")
	all.table <- all.table[order(all.table$cluster,adj.P.Val,-t,-logFC),]
	sig.table <- sig.table[order(sig.table$cluster,adj.P.Val,-t,-logFC),]
	if(!is.null(out.prefix))
	{
	    conn <- gzfile(sprintf("%s.limma.all.txt.gz",out.prefix),"w")
	    write.table(all.table,conn,row.names = F,quote = F,sep = "\t")
	    close(conn)
	    conn <- gzfile(sprintf("%s.limma.sig.txt.gz",out.prefix),"w")
	    write.table(sig.table,conn,row.names = F,quote = F,sep = "\t")
	    close(conn)
	}
    }

    return(list(all=all.table,sig=sig.table,fit=fit.list))
}







