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
#' @param assay.name assay name (default "exprs")
#' @details if x is an object of SingleCellExperiment, just clear the metadata;
#' if x is matrix/data.frame, convert it to an object of SingleCellExperiment.
#' Also a vector `display.name` can be provided, which would be used in some plots,
#' such as geneOnTSNE, heatmap. The row names of SingleCellExperiment object usually be gene id
#' (e.g. entrez ID, Ensemble ID), the `display.name` should be human readable gene name (
#' e.g. HGNC gene symbol). If `display.name` is NULL (default), the row names of SingleCellExperiment object
#' would be used.
#' @return an object of \code{SingleCellExperiment} class
#' @export
ssc.build <- function(x,display.name=NULL,assay.name="exprs")
{
  obj <- NULL
  if(class(x)=="SingleCellExperiment")
  {
    obj <- x
    metadata(obj)$ssc <- list()
  }else if(class(x) %in% c("matrix","data.frame")){
    obj <- SingleCellExperiment(assays = setNames(list(as.matrix(x)),assay.name))
  }else if(class(x) %in% c("dgCMatrix","dgTMatrix")){
    obj <- SingleCellExperiment(assays = setNames(list(x),assay.name))
  }
  metadata(obj)$ssc$colSet <- list()
  if(is.null(rowData(obj)[["display.name"]])){
    if(!is.null(display.name)){
      f.na <- is.na(display.name)
      display.name[f.na] <- row.names(obj)[f.na]
      rowData(obj)[,"display.name"] <- display.name
    }else{
      rowData(obj)[,"display.name"] <- row.names(obj)
    }
    if(is.null(names(rowData(obj)$display.name))){ names(rowData(obj)$display.name) <- row.names(obj) }
  }
  if(is.null(obj)){
    warning("SingleCellExperiment object building failed!")
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
    ret <- ret[display.name]
  }
  return(ret)
}

####### operations on sce object ######
#' sort by hierarchical clustering
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param order.col logical; wheter order columns (default: FALSE)
#' @param order.row logical; wheter order row (default: FALSE)
#' @param clustering.distance character; one of spearmn, pearson and euclidean (default: "spearman")
#' @param clustering.method character; method for hclust (default: "complete")
#' @param k.row integer; number of clusters in the rows (default: 1)
#' @param k.col integer; number of clusters in the columns (default: 1)
#' @importFrom stats hclust
#' @importFrom dendextend color_branches
#' @details order genes or cells according the assay, using hclust.
#' @export
ssc.assay.hclust <- function(obj,assay.name="exprs",
                             order.col=FALSE,order.row=FALSE,
                             clustering.distance="spearman",clustering.method="complete",
                             k.row=1,k.col=1)
{
    if(order.col && ncol(obj)>2)
    {
        branch.col <- FALSE
        obj.hclust.col <- NULL
        if(clustering.distance=="spearman" || clustering.distance=="pearson"){
            tryCatch({
                dist.out <- cor.BLAS(t(assay(obj,assay.name)),method=clustering.distance,nthreads=1)
                obj.hclust.col <- hclust(as.dist(1-dist.out), method=clustering.method)
            },error = function(e){
                cat("using spearman/pearson as distance failed;try to fall back to use euler distance ... \n");
            })
        }
        if(is.logical(branch.col) && !branch.col){
            obj.hclust.col <- hclust(dist(t(assay(obj,assay.name))),method=clustering.method)
        }
        branch.col <- dendextend::color_branches(as.dendrogram(obj.hclust.col),k=k.col)
        obj <- obj[,obj.hclust.col$order]
        metadata(obj)$assay.hclust$col <- obj.hclust.col
        metadata(obj)$assay.hclust$branch.col <- branch.col
    }
    if(order.row && nrow(obj)>2){
        branch.row <- FALSE
        obj.hclust.row <- NULL
        if(clustering.distance=="spearman" || clustering.distance=="pearson"){
            tryCatch({
                dist.out <- cor.BLAS(assay(obj,assay.name),method=clustering.distance,nthreads=1)
                obj.hclust.row <- hclust(as.dist(1-dist.out), method=clustering.method)
            },error = function(e){
                cat("using spearman/pearson as distance failed;try to fall back to use euler distance ... \n");
            })
        }
        if(is.logical(branch.row) && !branch.row){
            obj.hclust.row <- hclust(dist(assay(obj,assay.name)),method=clustering.method)
        }
        branch.row <- dendextend::color_branches(as.dendrogram(obj.hclust.row),k=k.row)
        obj <- obj[obj.hclust.row$order,]
        metadata(obj)$assay.hclust$row <- obj.hclust.row
        metadata(obj)$assay.hclust$branch.row <- branch.row
    }
    return(obj)
}

#' order genes and cells
#' @param obj object of \code{singleCellExperiment} class
#' @param columns.order character; columns of colData(obj) used for ordering (default: NULL)
#' @param gene.desc data.frame; it must contain columns geneID and Group (default: NULL)
#' @import data.table
#' @export
ssc.order <- function(obj,columns.order=NULL,gene.desc=NULL)
{
    if(!is.null(gene.desc) && ("Group" %in% colnames(gene.desc)) && ("geneID" %in% colnames(gene.desc))){
		gene.desc <- as.data.table(gene.desc)
		gene.desc <- gene.desc[order(Group),]
        obj <- obj[gene.desc$geneID,]
    }
    if(!is.null(columns.order)){
        annDF <- as.data.frame(colData(obj)[columns.order])
        annDF <- annDF[eval(parse(text=sprintf("with(annDF,order(%s))", columns.order))),,drop=F]
        obj <- obj[,rownames(annDF)]
    }
    return(obj)
}


#' calcualte the average expression of cells of each group
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param gene character; only consider the specified gnees (default: NULL)
#' @param column character; columns in colData(obj) to be averaged. (default: "majorCluster")
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param avg character; average method. can be one of "mean", "diff", "zscore" . (default: "mean")
#' @param ret.type character; return type. can be one of "data.melt", "data.cast", "data.mtx". (default: "data.melt")
#' @importFrom plyr ldply
#' @importFrom Matrix rowMeans
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowSds
#' @importFrom data.table dcast
#' @details multiple average methods are implemented
#' @export
ssc.average.cell <- function(obj,assay.name="exprs",gene=NULL,column="majorCluster",ncell.downsample=NULL,
                             avg="mean",ret.type="data.melt")
{
  if(!column %in% colnames(colData(obj))){
    warning(sprintf("column not in the obj: %s \n",column))
    return(NULL)
  }
  if(!is.null(gene)){
    obj <- obj[gene,]
  }

  #### downsample cells
  if(!is.null(ncell.downsample)){
	clust <- colData(obj)[,column]
	names(clust) <- colnames(obj)
	grp.list <- unique(clust)
	f.cell <- unlist(sapply(grp.list,function(x){
						 x <- names(clust[clust==x])
						 sample(x,min(length(x),ncell.downsample)) }))
	obj <- obj[,f.cell]
  }

  cls <- sort(unique(colData(obj)[,column]))
  data.melt.df <- ldply(cls,function(x){
    obj.in <- obj[,colData(obj)[,column]==x]
    avg.in <- NULL
    avg.in <- Matrix::rowMeans(assay(obj.in,assay.name))
    if(avg=="mean"){
      return(data.frame(geneID=names(avg.in),cls=x,avg=avg.in,
                        stringsAsFactors = F))
    }else if (avg=="diff"){
      obj.out <- obj[,colData(obj)[,column]!=x]
      avg.out <- Matrix::rowMeans(assay(obj.out,assay.name))
      return(data.frame(geneID=names(avg.out),cls=x,avg=avg.in-avg.out,
                        stringsAsFactors = F))
    }else if (avg=="zscore"){
      obj.out <- obj[,colData(obj)[,column]!=x]
      avg.out <- Matrix::rowMeans(assay(obj.out,assay.name))
      ##sd.r <- matrixStats::rowSds(assay(obj,assay.name))
      sd.r <- DelayedMatrixStats::rowSds(DelayedArray(assay(obj,assay.name)))
      dat.ret <- data.frame(geneID=names(avg.out),cls=x,avg=(avg.in-avg.out)/sd.r,
							stringsAsFactors = F)
	  dat.ret$avg[is.na(dat.ret$avg)] <- 0
	  return(dat.ret)
    }
  })
  if(ret.type=="data.melt"){
    return(data.melt.df)
  }else if(ret.type=="data.dcast"){
    dat.df <- dcast(data.melt.df,geneID~cls,value.var="avg")
    rownames(dat.df) <- dat.df[,1]
    dat.df <- dat.df[rownames(obj),]
    return(dat.df)
  }else if(ret.type=="data.mtx"){
    dat.df <- dcast(data.melt.df,geneID~cls,value.var="avg")
    dat.mtx <- as.matrix(dat.df[,-1])
    rownames(dat.mtx) <- dat.df[,1]
    dat.mtx <- dat.mtx[rownames(obj),]
    return(dat.mtx)
  }else if(ret.type=="sce"){
    dat.df <- dcast(data.melt.df,geneID~cls,value.var="avg")
    dat.mtx <- as.matrix(dat.df[,-1])
    rownames(dat.mtx) <- dat.df[,1]
    dat.mtx <- dat.mtx[rownames(obj),]
    obj.ret <- ssc.build(dat.mtx,assay.name=assay.name,display.name=rowData(obj)$display.name)
    if(is.factor(obj[[column]])){
        ##alevels <- intersect(levels(obj[[column]]),colnames(obj.ret))
        alevels <- levels(obj[[column]])
        colData(obj.ret)[,column] <- factor(colnames(obj.ret),levels=alevels)
    }else{
        colData(obj.ret)[,column] <- colnames(obj.ret)
    }
    return(obj.ret)
  }
}


#' scale the assay per gene
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param assay.new character; assay name to store the scaled- expression;if NULL, will be (assay.name).z (default: NULL)
#' @param covar character; perform scale in cells grouping by covar (default: "patient")
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
#' @param z.lo double; z-score lower boundary; if set, z-score lower than this will be set to this (default: NULL)
#' @param z.hi double; z-score higher boundary; if set, z-score higher than this will be set to this (default: NULL)
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table as.data.table
#' @importFrom plyr ldply
#' @export
ssc.assay.zscore <- function(obj,assay.name="exprs",assay.new=NULL,covar="patient",n.cores=NULL,
                             z.lo=NULL,z.hi=NULL)
{
    requireNamespace("data.table")

    dat.plot <- as.matrix(assay(obj,assay.name))

    cell.info <- as.data.table(colData(obj)[covar])
    cell.info$cellID <- colnames(obj)
    cell.info.split <- split(cell.info,by=covar,sorted=T)

    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = n.cores)
    dat.plot.z.df <- data.table(ldply(names(cell.info.split),function(x){
                dat.block <- dat.plot[,cell.info.split[[x]][["cellID"]],drop=F]
                #dat.block <- scale(t(dat.block))
                rowM <- rowMeans(dat.block, na.rm = T)
                rowSD <- apply(dat.block, 1, sd, na.rm = T)
                dat.block <- sweep(dat.block, 1, rowM)
                dat.block <- sweep(dat.block, 1, rowSD, "/")
                dat.block <- t(dat.block)
                dat.block.df <- cbind(cellID=rownames(dat.block),as.data.frame(dat.block),stringsAsFactors=F)
                return(dat.block.df)
            },.progress = "none",.parallel=T))
    dat.plot.z <- as.matrix(dat.plot.z.df[,-c("cellID")])
    rownames(dat.plot.z) <- dat.plot.z.df$cellID
    dat.plot.z <- t(dat.plot.z)

    if(!is.null(z.lo) && !is.null(z.hi)){
        dat.plot.z[dat.plot.z < z.lo] <- z.lo
        dat.plot.z[dat.plot.z > z.hi] <- z.hi
    }
    assay(obj,if(is.null(assay.new)) sprintf("%s.z",assay.name) else assay.new) <- dat.plot.z[,colnames(obj)]
    return(obj)
}

#######################################


#' Reduce dimension by various methods
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom stats cor prcomp
#' @importFrom BiocGenerics t
#' @importFrom uwot umap
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
                          reuse=F,ncore=NULL,seed=NULL,out.prefix=NULL)
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
                         out.prefix=out.prefix,n.cores=ncore,method=method.tsne)
        }
    }else if(method=="tsne"){
      proj_data <- run.tSNE(BiocGenerics::t(assay(obj[vgene,],assay.name)),tSNE.usePCA,tSNE.perplexity,
                            out.prefix=out.prefix,n.cores=ncore,method=method.tsne)
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
                   n.cores=ncore,method=method.tsne)
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
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param method.reduction character; which dimention reduction method to be used, should be one of
#' "iCor", "pca" and "none". (default: "iCor")
#' @param method character; clustering method to be used, should be one of "kmeans", "hclust", "dynamicTreeCut", "SNN", "dpclust", "adpclust" and "SC3". (default: "kmeans")
#' @param k.batch integer; number of clusters to be evaluated. (default: 2:6)
#' @param method.vgene character; variable gene identification method used. (default: "HVG.sd")
#' @param SNN.k integer; number of shared NN. (default: 10)
#' @param SNN.method character; cluster method applied on SNN， one of "greedy", "eigen", "infomap",
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
      clust.res <- leiden::leiden(igraph::as_adjacency_matrix(snn.gr),seed=myseed,...)
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
#' @importFrom SingleCellExperiment reducedDims
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

#' Plot on tSNE map
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param gene, character; genes to be showed. (default: NULL)
#' @param columns character; columns in colData(obj) to be showd. (default: NULL)
#' @param splitBy character; columns in colData(obj). Split the dataset to mupltiple subset then plot them one by one (default: NULL)
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
#' @param brewer.palette character; which palette to use. (default: "YlOrRd")
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param clamp integer vector; expression values will be clamped to the range defined by this parameter, such as c(0,15). (default: "none" )
#' @param do.scale logical; whether scale the expression value. (default: FALSE)
#' @param label double; label size. if NULL, no label showed. (default: NULL )
#' @param par.repel list; passed to geom_text_repel
#' @param par.geneOnTSNE character; other parameters of geneOnTSNE
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual theme_bw aes_string guides guide_legend coord_cartesian
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot save_plot plot_grid
#' @importFrom utils read.table
#' @importFrom RColorBrewer brewer.pal
#' @details If `gene` is not NULL, expression of the specified genes will be plot on the tSNE map; if columns in not
#' NULL, colData of obj with names in `columns` will be plot on the tSNE map. The tSNE map used is specified by option
#' `reduced.name` and `reduced.dim`. Both `gene` and `columns` can be non-NULL. For list `colSet`, each element define
#' a color mapping for the responding iterm in the `column`; if not specifed, automatically generated color mapping will
#' be used.
#' @export
ssc.plot.tsne <- function(obj, assay.name="exprs", gene=NULL, columns=NULL,splitBy=NULL,
                             plotDensity=F, colSet=list(),
                             reduced.name="iCor.tsne",reduced.dim=c(1,2),xlim=NULL,ylim=NULL,size=NULL,
                             brewer.palette="YlOrRd",adjB=NULL,clamp="none",do.scale=FALSE,
                             label=NULL,par.repel=list(force=1),
                             par.geneOnTSNE=list(scales="free",pt.order="value",pt.alpha=0.1),
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
            if(is.null(metadata(obj)$ssc$colSet[[cc]])){
                cc.values <- sort(unique(colData(obj)[,cc]))
                colSet[[cc]] <- structure(auto.colSet(length(cc.values),name = "Paired"),
                                          names=as.character(cc.values))
            }else{
                colSet[[cc]] <- metadata(obj)$ssc$colSet[[cc]]
            }
          }
          dat.plot <- data.frame(sample=rownames(dat.map),stringsAsFactors = F)
          if(!is.null(splitBy)){
            dat.plot <- as.data.frame(cbind(dat.plot,dat.map,colData(obj)[,c(cc,splitBy),drop=F]))
            colnames(dat.plot) <- c("sample","Dim1","Dim2",cc,"splitBy")
          }else{
            dat.plot <- as.data.frame(cbind(dat.plot,dat.map,colData(obj)[,cc,drop=F]))
            colnames(dat.plot) <- c("sample","Dim1","Dim2",cc)
          }
		  if(par.geneOnTSNE$pt.order=="value"){
			  dat.plot <- dat.plot[order(dat.plot[,cc]),]
		  }else if(par.geneOnTSNE$pt.order=="random"){
			  dat.plot <- dat.plot[sample(nrow(dat.plot),nrow(dat.plot)),]
		  }
          npts <- nrow(dat.plot)
          if(is.numeric(dat.plot[,cc])){
            nvalues <- Inf
            if(clamp!="none"){
                dat.plot[[cc]][ dat.plot[[cc]] < clamp[1] ] <- clamp[1]
                dat.plot[[cc]][ dat.plot[[cc]] > clamp[2] ] <- clamp[2]
            }
          }else{
            nvalues <- length(unique(dat.plot[,cc]))
          }
          p <- ggplot2::ggplot(dat.plot,aes(Dim1,Dim2)) +
            geom_point(aes_string(colour=cc),
                       show.legend=if(!is.numeric(dat.plot[,cc]) && nvalues>40) F else NA,
                       size=if(is.null(size)) auto.point.size(npts)*1.1 else size)
          if(!is.null(label)){
              dat.plot.label <- as.data.table(dat.plot)[,.(Dim1=median(.SD$Dim1),
                                                           Dim2=median(.SD$Dim2)),
                                                        by=cc]
              p <- p + do.call(ggrepel::geom_text_repel,c(list(aes_string("Dim1","Dim2",label = cc),
                                                               size=label,data=dat.plot.label),
                                                          par.repel))
          }
          if(!is.null(splitBy)){
            p <- p + ggplot2::facet_wrap(~splitBy)
          }
          if(is.numeric(dat.plot[,cc])){
            p <- p + scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, brewer.palette))
          }else{
            p <- p + scale_colour_manual(values = colSet[[cc]])
          }
          p <- p + theme_bw() + coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE) +
            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=if(nvalues<=26) 4 else 2.0),
                                                           label.theme = element_text(size=8)))

          return(p)
        })
        pp <- cowplot::plot_grid(plotlist=multi.p,ncol = if(length(columns)>1) 2 else 1,align = "hv")
        if(!is.null(out.prefix)){
          cowplot::save_plot(sprintf("%s.columnsOntSNE.pdf",out.prefix),pp,
                             ncol = if(length(columns)>1) 2 else 1,
                             base_aspect_ratio=base_aspect_ratio)
        }else{
          #print(pp)
          return(pp)
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

    dat.onTSNE <- assay(obj,assay.name)[gene,,drop=F]
    if(!is.null(adjB)){
      dat.onTSNE <- simple.removeBatchEffect(dat.onTSNE,batch=colData(obj)[[adjB]])
    }
    if(do.scale){
      dat.onTSNE <- t(scale(t(dat.onTSNE)))
    }
    p <- do.call(ggGeneOnTSNE,c(list(Y=dat.onTSNE, dat.map=dat.map, gene.to.show=gene,
                                     p.ncol=p.ncol,xlim=xlim,ylim=ylim,
                                     size=size,width=width,height=height,
                                     clamp=clamp,
                                     out.prefix=out.prefix),
                                par.geneOnTSNE))
    if(is.null(out.prefix)){
      #print(p)
      return(p)
    }
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
#' @param clamp integer vector; expression values will be clamped to the range defined by this parameter. (default: c(0,12))
#' @param out.prefix character; output prefix. (default: NULL)
#' @param p.ncol integer; number of columns in the figure layout. (default: 3)
#' @param base_aspect_ratio numeric; base_aspect_ratio, used for plotting metadata. (default 1.1)
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param do.scale logical; whether scale the expression value. (default: FALSE)
#' @param ... parameter passed to cowplot::save_plot
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_violin scale_fill_gradient2 theme_bw theme aes_string facet_grid element_text
#' @importFrom cowplot save_plot plot_grid
#' @importFrom data.table melt
#' @details If `gene` is not NULL, violin of the genes' expression will be plot; if columns in not
#' NULL, colData of obj with names in `columns` will be plot in violin.
#' @export
ssc.plot.violin <- function(obj, assay.name="exprs", gene=NULL, columns=NULL,
                            group.var="majorCluster",clamp=c(0,12),adjB=NULL,do.scale=F,
                            out.prefix=NULL,p.ncol=1,base_aspect_ratio=1.1,...)
{
  requireNamespace("ggplot2")
  requireNamespace("data.table")
  if(!is.null(gene))
  {
	  gene <- ssc.displayName2id(obj,display.name = gene)

	  dat.violin <- assay(obj,assay.name)[gene,,drop=F]
	  if(!is.null(adjB)){
		dat.violin <- simple.removeBatchEffect(dat.violin,batch=colData(obj)[[adjB]])
	  }
      if(do.scale){
        dat.violin <- t(scale(t(dat.violin)))
      }

	  dat.plot <- as.matrix(t(dat.violin))
	  colnames(dat.plot) <- ssc.id2displayName(obj,colnames(dat.plot))
	  dat.plot.df <- data.table::data.table(sample=rownames(dat.plot),stringsAsFactors = F)
	  dat.plot.df <- cbind(dat.plot.df,as.data.frame(colData(obj)[,group.var,drop=F]))
	  dat.plot.df <- cbind(dat.plot.df,dat.plot)
	  dat.plot.df <- data.table::melt(dat.plot.df,id.vars=c("sample",group.var),
									  variable.name="gene",value.name=assay.name)
	  dat.plot.df.grpMean <- dat.plot.df[,lapply(.SD,mean),by=c("gene",group.var),.SDcols=assay.name]
	  colnames(dat.plot.df.grpMean) <- c("gene",group.var,"meanExp")
	  dat.plot.df <- dat.plot.df.grpMean[dat.plot.df,,on=c("gene",group.var)]
      dat.plot.df[,gene:=factor(gene,levels=colnames(dat.plot),ordered=T)]

      if(is.null(clamp)){
          clamp <- quantile(dat.plot.df[[assay.name]],c(0.05,0.95))
      }
      dat.plot.df[meanExp<clamp[1],meanExp:=clamp[1],]
      dat.plot.df[meanExp>clamp[2],meanExp:=clamp[2],]
      dat.plot.df[[assay.name]][ dat.plot.df[[assay.name]] < clamp[1] ] <- clamp[1]
      dat.plot.df[[assay.name]][ dat.plot.df[[assay.name]] > clamp[2] ] <- clamp[2]

	  p <- ggplot(dat.plot.df, aes_string(group.var[1], assay.name))
	  if(length(group.var)==1){
		p <- p +
		  geom_violin(scale = "width",aes(fill=meanExp),color=NA,show.legend = T) +
		  scale_fill_gradient2(low = "yellow",mid = "red",high = "black",midpoint = mean(clamp), limits=clamp)
	  }else if(length(group.var)==2)
	  {
		p <- p +
			geom_boxplot(aes_string(colour = group.var[2])) +
			scale_colour_brewer(palette = "Set1")
			#geom_violin(scale = "width",aes_string(fill="meanExp",linetype=group.var[2],color=group.var[2]),
			#            show.legend = T) +
			#scale_fill_gradient2(low = "yellow",mid = "red",high = "black",midpoint = mean(clamp),limits=clamp)
	  }
	  p <- p +
		theme_bw(base_size = 12) +
#		facet_grid(gene ~ .,switch = "y",scales = "free_y") +
        facet_wrap(gene~.,strip.position = "left",scales="free_y",dir="v",ncol=p.ncol) +
		theme(axis.text.x = element_text(angle = 60, hjust = 1),
              strip.placement = "inside")
  } else if(!is.null(columns)){
	  dat.plot.df <- as.data.table(cbind(data.frame(cellID=colnames(obj),stringsAsFactors=F),
						  as.data.frame(colData(obj)[,c(group.var,columns),drop=F])))
	  dat.plot.df <- melt(dat.plot.df,id.vars=c("cellID",group.var))
	  p <- ggplot(dat.plot.df, aes_string(group.var[1], "value"))
      if(length(group.var)==1){
          p <- p + geom_boxplot()
      }else if(length(group.var)==2){
          p <- p + geom_boxplot(aes_string(colour=group.var[2]))
          if(is.null(metadata(obj)$ssc$colSet)){
              p <- p + scale_colour_brewer(palette = "Set1")
          }else{
              p <- p + scale_colour_manual(values = metadata(obj)$ssc$colSet[[group.var[2]]])
          }
      }
      p <- p + theme_bw(base_size = 12) +
			facet_grid(variable ~ .,switch = "y",scales = "free_y") +
			theme(axis.text.x = element_text(angle = 60, hjust = 1),strip.placement = "inside")
  }
  if(!is.null(out.prefix)){
    cowplot::save_plot(sprintf("%s.violin.%s.pdf",out.prefix,if(!is.null(gene)) "gene" else "columns"),p,
                       ncol = p.ncol,
                       base_aspect_ratio=base_aspect_ratio,...)
  }else{
    return(p)
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
#' @param verbose logical; whether output all genes' result. (default: FALSE)
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply
#' @importFrom dplyr inner_join
#' @details identify marker genes based on aov and AUC.
#' @export
ssc.clusterMarkerGene <- function(obj, assay.name="exprs", ncell.downsample=NULL,
                                  group.var="majorCluster",batch=NULL,
                                  assay.bin=NULL, out.prefix=NULL,n.cores=NULL, do.plot=T,
                                  F.FDR.THRESHOLD=0.01,pairwise.P.THRESHOLD=0.01,pairwise.FC.THRESHOLD=1,
                                  verbose=F)
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
    if(length(unique(clust))<2 || is.null(obj) || !all(table(clust) > 1)){
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
        getAUC(dat.to.test[x,],clust,use.rank = F,geneID=x)
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
                                 out.prefix=out.prefix,mod=NULL,
                                 F.FDR.THRESHOLD=F.FDR.THRESHOLD,
                                 HSD.FDR.THRESHOLD=pairwise.P.THRESHOLD,
                                 HSD.FC.THRESHOLD=pairwise.FC.THRESHOLD,
                                 verbose=verbose,n.cores=n.cores, gid.mapping=gid.mapping)
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
    f.isSig <- intersect(rownames(.aov.res$aov.out.sig),rownames(.gene.table))
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
#' @param verbose integer; verbose level. (default: 0)
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom limma lmFit eBayes topTable
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply llply
#' @details identify differential genes using limma
#' @export
ssc.DEGene.limma <- function(obj, assay.name="exprs", ncell.downsample=NULL,
                                  group.var="majorCluster",group.list=NULL,group.mode="multi",batch=NULL,
                                  out.prefix=NULL,n.cores=NULL, do.plot=T,
                                  T.fdr=0.01,T.logFC=1,
                                  verbose=0)
{
#    requireNamespace("doParallel")
    requireNamespace("plyr")
    requireNamespace("dplyr")

    if(!is.null(ncell.downsample)){
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
    if(!all(stat.clust[group.list]>=2)){
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
		if(group.mode=="multiAsTwo"){
			xlabel[xlabel!=x] <- "_control"
			xlabel[xlabel==x] <- "_case"
			xgroup <- sprintf("_case:%s",x)
		}
        out.limma <- run.limma.matrix(assay(obj,assay.name),xlabel,batch=batchV,
                                      out.prefix=if(is.null(out.prefix)) NULL else sprintf("%s.%s",out.prefix,x),
                                      group=xgroup,
                                      T.fdr=T.fdr,T.logFC=T.logFC,verbose=verbose,n.cores=1,
                                      gid.mapping=gid.mapping, do.voom=F)

    },.parallel=T)
    names(out) <- group.list

    all.table <- data.table(ldply(group.list,function(x){ out[[as.character(x)]]$all }))
    sig.table <- data.table(ldply(group.list,function(x){ out[[as.character(x)]]$sig }))

	if(verbose>2){
		res.aov <- findDEGenesByAOV(assay(obj,assay.name),clust,batch=batchV, out.prefix=NULL,
									n.cores=n.cores, gid.mapping=gid.mapping)
		res.aov$aov.out$F.rank <- rank(-res.aov$aov.out$F)/nrow(res.aov$aov.out)
		all.table <- merge(all.table,res.aov$aov.out[,c("geneID","F","F.pvalue","F.adjp","F.rank")],by="geneID")
		sig.table <- merge(sig.table,res.aov$aov.out[,c("geneID","F","F.pvalue","F.adjp","F.rank")],by="geneID")
		all.table <- all.table[order(cluster,adj.P.Val,-t,-logFC),]
		sig.table <- sig.table[order(cluster,adj.P.Val,-t,-logFC),]
	    if(!is.null(out.prefix))
		{
			write.table(all.table,sprintf("%s.limma.all.txt",out.prefix),row.names = F,quote = F,sep = "\t")
			write.table(sig.table,sprintf("%s.limma.sig.txt",out.prefix),row.names = F,quote = F,sep = "\t")
		}
	}

    return(list(all=all.table,sig=sig.table))
}


#' plot heatmap
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param out.prefix character; output prefix. (default: NULL)
#' @param ncell.downsample integer; number of cells downsample to. (default: NULL)
#' @param ave.by character; average the expression profile grouping by this. (default: NULL)
#' @param columns character; columns in colData(obj) to be showd. must be subset of columns of colData(obj) and ave.by (if it's not NULL) (default: NULL)
#' @param columns.order character; columns of colData(obj) used for ordering (default: NULL)
#' @param gene.desc data.frame; it must contain columns geneID and Group (default: NULL)
#' @param colSet list; mapping iterms in the names to colors in the values. (default: list())
#' @param pdf.width double; width of the pdf file. (default:16)
#' @param pdf.height double; height of the pdf file. (default:15)
#' @param do.scale logical; wheter scale the rows, just for visualization. (default: TRUE)
#' @param z.lo double; z-score lower boundary; z-score lower than this will be set to this (default: -2.5)
#' @param z.hi double; z-score higher boundary; z-score higher than this will be set to this (default: 2.5)
#' @param z.step double; z-score step, used for coloring the expression value (default: 1)
#' @param exp.title character; title for the expression legend (default: "Exp")
#' @param do.clustering.row logical; wheter order row (default: TRUE)
#' @param do.clustering.col logical; wheter order columns (default: TRUE)
#' @param dend.col dendrogram of the columns, 'cluster_columns' of ComplexHeatmap::Heatmap  (default: FALSE)
#' @param dend.row dendrogram of the rows, 'cluster_rows' of ComplexHeatmap::Heatmap (default: FALSE)
#' @param clustering.distance character; one of spearmn, pearson and euclidean (default: "spearman")
#' @param clustering.method character; method for hclust (default: "complete")
#' @param k.row integer; number of clusters in the rows (default: 1)
#' @param k.col integer; number of clusters in the columns (default: 1)
#' @param palette.name character; which palette to use, such as "RdBu","RdYlBu" (default: NULL)
#' @param row.split vector; used for row; (default: NULL)
#' @param column.split vector; used for column; (default: NULL)
#' @param annotation_legend_param list; (default: list())
#' @param ann.bar.height double; height of the top annotation (default: 1.5)
#' @param mytitle character; title of the figure (default: "")
#' @param ... parameters pass to sscClust:::plot.matrix.simple()
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap decorate_annotation
#' @importFrom circlize colorRamp2
#' @importFrom gridBase baseViewports
#' @importFrom grid pushViewport grid.text gpar unit
#' @importFrom RColorBrewer brewer.pal
#' @details identify marker genes based on aov and AUC.
#' @export
ssc.plot.heatmap <- function(obj, assay.name="exprs",out.prefix=NULL,
                             ncell.downsample=NULL,ave.by=NULL,
                             columns=NULL,columns.order=NULL,gene.desc=NULL,
                             colSet=list(), pdf.width=16,pdf.height=15,
                             do.scale=TRUE,z.lo=-2.5,z.hi=2.5,z.step=1,exp.title="Exp",
                             do.clustering.row=T,do.clustering.col=T,
                             dend.col=FALSE,dend.row=FALSE,
                             clustering.distance="spearman",clustering.method="complete",k.row=1,k.col=1,
                             palette.name=NULL,row.split=NULL,column.split=NULL,
                             annotation_legend_param=list(),ann.bar.height=1.5, mytitle="",...)
{
    requireNamespace("ComplexHeatmap")
    requireNamespace("circlize")
    requireNamespace("gridBase")
    requireNamespace("grid")
    requireNamespace("RColorBrewer")

    if(!is.null(gene.desc) && ("Group" %in% colnames(gene.desc)) && ("geneID" %in% colnames(gene.desc))){
        obj <- obj[gene.desc$geneID,]
    }
    if(!is.null(ncell.downsample) && ncell.downsample < ncol(obj) ){
        obj <- obj[,sample(seq_len(ncol(obj)),ncell.downsample)]
    }

    n <- nrow(obj)
    m <- ncol(obj)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }

    if (!is.null(row.split) && is.null(names(row.split))) {
        names(row.split) <- rownames(obj)
    }
    if (!is.null(column.split) && is.null(names(column.split))) {
        names(column.split) <- colnames(obj)
    }

    #### sort
    if(is.null(ave.by)){
        obj <- ssc.assay.hclust(obj,assay.name=assay.name,
                                order.col=if(is.logical(dend.col) && FALSE==dend.col) do.clustering.col else FALSE,
                                order.row=if(is.logical(dend.row) && FALSE==dend.row) do.clustering.row else FALSE,
                                clustering.distance="spearman",clustering.method="complete",
                                k.row=1,k.col=1)
    }else{
        avg.colDat <- unique((colData(obj)[,unique(c(ave.by,columns,columns.order)),drop=F]))
        obj <- ssc.average.cell(obj,assay.name=assay.name,column=ave.by,ret.type="sce")
        #columns <- intersect(ave.by,columns)
        #columns.order <- intersect(ave.by,columns.order)
        if(nrow(avg.colDat)==ncol(obj) && all(avg.colDat[[ave.by]] %in% colnames(obj)) ){
            rownames(avg.colDat) <- avg.colDat[[ave.by]]
            colData(obj) <- avg.colDat[colnames(obj),,drop=F]
        }
        obj <- ssc.assay.hclust(obj,assay.name,order.col=do.clustering.col,order.row=do.clustering.row)
    }

    #### visualization of annotation on top of heatmap
    ha.col <- NULL
    annDF <- data.frame()
    if(!is.null(columns))
    {
        if(!is.null(columns.order)){
            obj <- ssc.order(obj,columns.order=columns.order)
        }
        annDF <- as.data.frame(colData(obj)[columns])
        if(length(colSet)==0) {
            for(i in seq_along(columns)){
                x <- columns[i]
                if(class(colData(obj)[,x])=="numeric"){
                    if(all(colData(obj)[,x]<=1) && all(colData(obj)[,x]>=0)){
                        Y.level <- c(0,1)
                    }else{
                        Y.level <- pretty(colData(obj)[,x],n=8)
                    }
                    # continious version
                    colSet[[x]] <- colorRamp2(seq(Y.level[1],Y.level[length(Y.level)],length=7),
                                              rev(brewer.pal(n = 7, name = "RdYlBu")),
                                              space="LAB")
                    annotation_legend_param[[x]] <- list(color_bar="continuous",
                                                         legend_direction="horizontal",
                                                         legend_width=unit(4, "cm"),
                                                         legend_height=unit(2, "cm"))
                }else{
                    group.value <- sort(unique(colData(obj)[,x]))
                    colSet[[x]] <- structure(sscClust:::auto.colSet(length(group.value),name="Accent"),
                                             names=as.character(group.value))
                }
            }
        }

        g.show.legend <- T
        ha.col <- ComplexHeatmap::HeatmapAnnotation(df = annDF, col = colSet,
                                    show_legend = g.show.legend,
									simple_anno_size = unit(ann.bar.height, "cm"),
                                    annotation_legend_param = annotation_legend_param)
        ###top_annotation_height <- unit(ann.bar.height * ncol(annDF), "cm")
    }

    obj <- ssc.order(obj,columns.order=NULL,gene.desc=gene.desc)

    #### scale data for visualization
    dat.plot <- as.matrix(assay(obj,assay.name))
    rownames(dat.plot) <- unname(rowData(obj)$display.name)
    #### scale by row
    if(do.scale)
    {
        rowM <- rowMeans(dat.plot, na.rm = T)
        rowSD <- apply(dat.plot, 1, sd, na.rm = T)
        dat.plot <- sweep(dat.plot, 1, rowM)
        dat.plot <- sweep(dat.plot, 1, rowSD, "/")
        if(!is.null(z.lo)){ dat.plot[dat.plot < z.lo] <- z.lo }
        if(!is.null(z.hi)){ dat.plot[dat.plot > z.hi] <- z.hi }
    }else{
        ###tmp.var <- pretty(abs(dat.plot),n=8)
        tmp.var <- pretty((dat.plot),n=8)
        if(is.null(z.lo)){ z.lo <- tmp.var[1] }
        if(is.null(z.hi)) { z.hi <- tmp.var[length(tmp.var)] }
        if(is.null(z.step)) { z.step <- tmp.var[2]-tmp.var[1] }
    }

    ##### plot
    if(!is.null(out.prefix)){ pdf(sprintf("%s.pdf",out.prefix),width=pdf.width,height=pdf.height) }
    par(mar=c(4,12,4,4))
    plot.new()
    title(main = mytitle,cex.main=2)
    #legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)

    ### Integrating Grid Graphics Output with Base Graphics Output
    vps <- gridBase::baseViewports()
    grid::pushViewport(vps$inner, vps$figure, vps$plot)

    if(is.null(palette.name)){
        exp.palette <- rev(brewer.pal(n = 7, name = ifelse(do.scale,"RdBu","RdYlBu")))
    }else{
        exp.palette <- rev(brewer.pal(n = 7, name = palette.name))
    }

    if(!is.null(row.split)){
        row.split <- row.split[rownames(dat.plot)]
    }
    if(!is.null(column.split)){
        column.split <- column.split[colnames(dat.plot)]
    }
    ht <- plot.matrix.simple(dat.plot,out.prefix=NULL,exp.name=exp.title,show.number=F,
                               do.clust=NULL,z.lo=z.lo,z.hi=z.hi,palatte=exp.palette,
                               clust.row=FALSE,clust.column=FALSE,show.dendrogram=FALSE,
                               returnHT=TRUE,column_split=column.split,
                               par.legend=list(at = seq(z.lo,z.hi,z.step)),
                               top_annotation = ha.col,...)

    ComplexHeatmap::draw(ht, newpage= FALSE,merge_legends = TRUE,split=row.split)
    if(!is.null(ha.col)){
        for(i in seq_along(names(ha.col@anno_list))){
          ComplexHeatmap::decorate_annotation(names(ha.col@anno_list)[i],
                                {grid.text(names(ha.col@anno_list)[i], unit(-4, "mm"),
                                           gp=grid::gpar(fontsize=14),just = "right")})
        }
    }
    if(!is.null(out.prefix)){ dev.off() }
}

#' plot gene expression density
#' @param obj object of \code{singleCellExperiment} class
#' @param out.prefix character; output prefix. required
#' @param gene.id should be in rownames(obj); genes to plot
#' @param gene.symbol should be in rowData(obj)[,"display.name"]; genes to plot
#' @param assay.name character; which assay (default: "exprs")
#' @param pallete.name character; pallete to use. (default: "heat")
#' @param expT double; expression threshold of the genes. (default: c(0.3,0.3))
#' @param ann.txt.dis double; adjust the position of the annotation text (default: 0.3)
#' @param ann.txt.cex double; cex for annotation text (default: 1.2)
#' @param my.title character; title of the figure (default: "")
#' @importFrom ks kde
#' @importFrom fields image.plot
#' @importFrom scales viridis_pal
#' @details make density plot of genes. Note, density estimation from ggplot2 is different (and not so 'effective' as that from ks::kde). One of gene.id and gene.symbol must be provided.
#' @export
ssc.plotGeneDensity <- function(obj,out.prefix,gene.id,gene.symbol,assay.name="norm_exprs",
								expT=c(0.3,0.3),pallete.name="heat",
								#adjB=NULL,do.scale=F,
								ann.txt.dis=0.3,ann.txt.cex=1.2,my.title="")
{

	if(missing(gene.id) && missing(gene.symbol)){
		warning("No gene.id or gene.symbol provided!")
		return(NULL)
	}
	if(missing(gene.id) && !is.null(gene.symbol)){
		gene.list <- rowData(obj)[,"display.name"][which(rowData(obj)[,"display.name"] %in% gene.symbol)]
	}else{
		gene.id <- intersect(rownames(obj),gene.id)
		gene.list <- rowData(obj)[,"display.name"][gene.id]
	}
	if(length(gene.list)==0){ 
		warning("No data found for the provided genes!")
		return(NULL)
	}
	obj <- obj[names(gene.list),]
	dat.block <- as.matrix(assay(obj,sprintf("%s",assay.name)))
	dat.plot <- cbind(data.table(cellID=colnames(obj)),(t(dat.block)))
	gene.list <- rowData(obj)[,"display.name"]

	if(length(gene.list)==2){
		colnames(dat.plot)[c(2,3)] <- c("x","y")
		par.old = par(no.readonly = T)

		.density <- ks::kde(dat.plot[,c("x","y")])
		.density.x <- density(dat.plot$x)
		.density.y <- density(dat.plot$y)
		zz <- c(5,10,20,30,40,50,60,70,80,90,95)
		if(pallete.name=="heat"){
			col.den2d <- c("transparent", rev(heat.colors(length(zz))))
		}else if(pallete.name=="rainbow"){
			col.den2d <- c("transparent", rev(rainbow(length(zz), end = 4/6)))
		}else if(pallete.name=="viridis"){
			col.den2d <- c("transparent", (scales::viridis_pal()(length(zz))))
		}

		dat.plot.x.range <- pretty(range(dat.plot$x))
		dat.plot.y.range <- pretty(range(dat.plot$y))
		dat.plot.x.range <- c(dat.plot.x.range[1],dat.plot.x.range[length(dat.plot.x.range)])
		dat.plot.y.range <- c(dat.plot.y.range[1],dat.plot.y.range[length(dat.plot.y.range)])

		mar.left <- 7
		pdf(sprintf("%s.dens2D.pdf",out.prefix),width = 8,height = 7)
		par(mar = c(0.8,mar.left,1,0),cex.lab=2,cex.axis=1.5)
		#layout(matrix(1:6, nrow = 2, byrow = T), widths = c(10,3,2.5), heights = c(3,10))
		nf <- layout(matrix(c(4,4,4,1,0,0,2,3,5),3,3,byrow = TRUE), c(10,3,2.5), c(1.5,3,10))
		#layout.show(nf)
		#dev.off()

		### upper
		plot(NULL, type = "n",  ylab = "",xlab="", xlim = dat.plot.x.range, ylim = c(0, max(.density.x$y)),
			 main = NA, axes = F, xaxs="i")
		lines(.density.x, col = "darkblue", lwd = 2)
		abline(v=expT[1],col="red4",lty=2,lwd=1.5)
		title(ylab = "density", line = 4.5)
		axis(2, las = 1)

		### main
		par(mar = c(6,mar.left,0,0))
	#	plot(NULL, type = "n",  ylab = "", xlim = dat.plot.x.range, ylim = dat.plot.y.range,
	#		 main = NA, axes = F, xaxs="i",yaxs="i",xlab="")
		plot(.density,display="filled.contour2", cont=zz,xlab="", ylab="",col=col.den2d,
			 xlim=dat.plot.x.range,ylim=dat.plot.y.range)
		title(ylab = gene.list[2], line = 4.5)
		title(xlab = gene.list[1])
		.addAnn <- function()
		{
			abline(v=expT[1],col="red4",lty=2,lwd=1.5)
			abline(h=expT[2],col="red4",lty=2,lwd=1.5)
			nn <- sum(dat.plot$x >= expT[1] & dat.plot$y >= expT[2])
			##mtext(text = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),side = 3,line=-2.0,adj = 0.95,xpd=T)
			text(par('usr')[2],par('usr')[4]-ann.txt.dis,
				 labels = sprintf("%4.2f %%\n(%d/%d)", nn*100/nrow(dat.plot),nn,nrow(dat.plot)),
				 adj = 1.1,xpd=T,cex=ann.txt.cex)
			nn <- sum(dat.plot$x < expT[1] & dat.plot$y >= expT[2])
			##mtext(text = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),side = 3,line=-2.0,adj = 0.05,xpd=T)
			text(expT[1],par('usr')[4]-ann.txt.dis,
				 labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),
				 adj = 1.1,xpd=T,cex=ann.txt.cex)
			nn <- sum(dat.plot$x < expT[1] & dat.plot$y < expT[2])
			text(expT[1],expT[2]-ann.txt.dis,
				 labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),
				 adj=1.1,xpd=T,cex=ann.txt.cex)
			nn <- sum(dat.plot$x >= expT[1] & dat.plot$y < expT[2])
			text(par('usr')[2],expT[2]-ann.txt.dis,
				 labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),
				 adj=1.1,xpd=T,cex=ann.txt.cex)
	#        mtext(text=sprintf("Pearson Cor: %4.2f (p value: %4.2e)",
	#						   cor.pcc.out$estimate,cor.pcc.out$p.value),
	#			  side=3,line=2.5,adj=0.5)
	#        mtext(text=sprintf("Spearman Cor: %4.2f (p value: %4.2e)",
	#						   cor.spe.out$estimate,cor.spe.out$p.value),
	#			  side=3,line=1.5,adj=0.5)
	#        mtext(text=sprintf("%s",sample.id),side=3,line=0.5,adj=0.5)
		}
		.addAnn()
		box()

		# right density plot
		par(mar = c(6,0.8,0,1))
		plot(NULL, type = "n", xlab = "density",ylab="",
			 ylim = dat.plot.y.range, xlim = c(0, max(.density.y$y)),
			 main = NA, axes = F, yaxs="i")
		lines(x=.density.y$y,y = .density.y$x, col = "darkblue", lwd = 2)
		abline(h=expT[2],col="red4",lty=2,lwd=1.5)
		axis(1,hadj = 0)

		### title
		par(mar=c(0,mar.left,0,4))
		plot.new()
		mtext(text=my.title,side=3,cex=1.8,line=-3,adj=0.5)
		
		### legend
		par(mar = c(6,0.25,0,1))
		plot.new()
		fields::image.plot(zlim=c(0,zz[length(zz)]),legend.only=TRUE,add=T,
				   col = col.den2d,
				   axis.args=list( at=zz, labels=sprintf("%s%%",100-zz)),
				   legend.width=12,legend.mar=28)
		par(par.old)
		dev.off()

		###### density estimation from ggplot2 is different (and not so 'effective' as that from ks::kde)
		#	p <- ggplot(dat.plot, aes(x, y)) +
		#		geom_point(alpha=0.5,size=0,color="white") +
		#		#stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE,bins=150) +
		#		#geom_density2d(bins=150) +
		#		#ggalt::stat_bkde2d(bandwidth=c(0.1,0.1),aes(fill = ..nlevel..), geom = "polygon")+
		#        scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))+
		#		#geom_contour() +
		#		#guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10)) +
		#		theme_bw()
		#	p2 <- ggExtra::ggMarginal(p, type = 'density')
		#	ggsave(sprintf("%s.test.00.pdf",out.prefix),p2,width=6,height=4)
	}
}

#' downsample
#' @param obj object of \code{singleCellExperiment} class
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param group.var character; column in the colData(obj) used for grouping. (default: "majorCluster")
#' @param rn.seed integer; random number seed (default: 9999)
#' @details return a downsampled object of \code{singleCellExperiment} class.
#' @export
ssc.downsample <- function(obj, ncell.downsample=NULL, group.var="majorCluster",rn.seed=9999)
{
    #### downsample cells
    set.seed(rn.seed)
    if(!is.null(ncell.downsample)){
        clust <- colData(obj)[,group.var]
        names(clust) <- colnames(obj)
        grp.list <- unique(clust)
        f.cell <- unlist(sapply(grp.list,function(x){
                             x <- names(clust[clust==x])
                             sample(x,min(length(x),ncell.downsample)) }))
        obj <- obj[,f.cell]
    }
    return(obj)
}


#' Add module scores for gene expression programs in single cells. Mofidy from Seurat::AddModuleScore
#'
#' Calculate the average expression levels of each program (cluster) on single cell level,
#' subtracted by the aggregated expression of control feature sets.
#' All analyzed features are binned based on averaged expression, and the control features are
#' randomly selected from each bin.
#'
#' @param obj object of SingleCellExperiment
#' @param features Feature expression programs in named list
#' @param pool List of features to check expression levels agains, defaults to \code{rownames(x = object)}
#' @param nbin Number of bins of aggregate expression levels for all analyzed features
#' @param ctrl Number of control features selected from the same bin per analyzed feature
#' @param assay.name Name of assay to use
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param do.scale logical; scale the data (default: TRUE)
#' @param seed Set a random seed
#' @return Returns a SingleCellExperiment object with module scores added to object meta data
#' @importFrom ggplot2 cut_number
#' @importFrom Matrix rowMeans colMeans
#' @references Tirosh et al, Science (2016); Seurat's source code
#' @export
ssc.moduleScore <- function(obj, features, pool = NULL,
							nbin = 24, ctrl = 100, assay.name = "exprs",
							adjB=NULL,do.scale=T,seed = 1)
{
	set.seed(seed = seed)
	if(is.null(names(features))){
		names(features) <- sprintf("M%02d",seq_along(features))
	}
	features <- llply(features,function(x){ intersect(x,rowData(obj)$display.name) })
	if(any(sapply(features,length)<1)){
		warning(sprintf("no expression data of the provided genes!\n"))
		return(obj)
	}
	features <- llply(features,function(x){ rownames(obj)[match(x,rowData(obj)$display.name)] })
	if(is.null(pool)){
		pool <- rownames(obj)
	}else{
		pool <- rownames(obj)[match(pool,rowData(obj)$display.name)]
	}
	assay.data <- assay(obj,assay.name)
	data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
	data.avg <- data.avg[order(data.avg)]
	data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30,
						   n = nbin, labels = FALSE, right = FALSE)
	###data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg) / (nbin + 1))))
	names(x = data.cut) <- names(x = data.avg)

	ctrl.use <- llply(features,function(x){
							  feat.rnd <- c()
							  for(i in seq_along(x)){
									feat.rnd <- c(feat.rnd,
												  names(sample(data.cut[which(data.cut==data.cut[x[i]])],size=ctrl,replace=F))
												  )
							  }
							  return(feat.rnd)
						   })
    .calScore <- function(x){
		 #### adjB and  scale
		 if(!is.null(adjB)){
			 dat.use <- simple.removeBatchEffect(assay.data[x,,drop=F],batch = obj[[adjB]])
		 }else{
			 dat.use <- t(scale(t(assay.data[x,,drop=F]),scale=F))
		 }
		 if(do.scale){
			 dat.use <- t(scale(t(dat.use)))
		 }
		 ###
		 Matrix::colMeans(dat.use)
	}
	ctrl.scores <- laply(ctrl.use,.calScore,.drop=F)
	features.scores <- laply(features,.calScore,.drop=F)
	features.scores.use <- features.scores - ctrl.scores
	rownames(features.scores.use) <- names(features)
	features.scores.use <- as.data.frame(x = t(x = features.scores.use))
	obj[[colnames(x = features.scores.use)]] <- features.scores.use[[1]]
	return(obj)
}

#' scale the data
#' @param obj object of \code{singleCellExperiment} class
#' @param gene.id should be in rownames(obj); genes to plot
#' @param gene.symbol should be in rowData(obj)[,"display.name"]; genes to plot
#' @param assay.name character; which assay (default: "exprs")
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param do.scale logical; whether scale the data. (default: F)
#' @details scale the data specified in assay.name, then store the scaled data to ${assay.name}.scale. One of gene.id and gene.symbol must be provided.
#' @export
ssc.scale <- function(obj,gene.id,gene.symbol,assay.name="norm_exprs",adjB=NULL,do.scale=F)
{
	if(missing(gene.id) && missing(gene.symbol)){
		warning("No gene.id or gene.symbol provided!")
		return(NULL)
	}
	if(missing(gene.id) && !is.null(gene.symbol)){
		gene.list <- rowData(obj)[,"display.name"][which(rowData(obj)[,"display.name"] %in% gene.symbol)]
	}else{
		gene.id <- intersect(rownames(obj),gene.id)
		gene.list <- rowData(obj)[,"display.name"][gene.id]
	}
	if(length(gene.list)==0){ 
		warning("No data found for the provided genes!")
		return(NULL)
	}
	obj <- obj[names(gene.list),]
	dat.block <- assay(obj,assay.name)
	if(!is.null(adjB)){
		dat.block <- simple.removeBatchEffect(dat.block,batch=obj[[adjB]])
	}
	if(do.scale){
		dat.block <- t(scale(t(dat.block)))
	}
	assay(obj,sprintf("%s.scale",assay.name)) <- dat.block
	return(obj)
}

