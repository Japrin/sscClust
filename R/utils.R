#' Correlation calculation. use BLAS and data.table to speed up.
#'
#' @importFrom data.table frank
#' @importFrom RhpcBLASctl omp_get_num_procs omp_set_num_threads
#' @param x matrix; input data, rows for variable (genes), columns for observations (cells).
#' @param method character; method used. (default: "pearson")
#' @param nthreads integer; number of threads to use. if NULL, automatically detect the number. (default: NULL)
#' @details calcualte the correlation among variables(rows)
#' @return correlation coefficient matrix among rows
cor.BLAS <- function(x,method="pearson",nthreads=NULL)
{
  if(is.null(nthreads))
  {
    nprocs <- RhpcBLASctl::omp_get_num_procs()
    RhpcBLASctl::omp_set_num_threads(max(nprocs-1,1))
  }else{
    RhpcBLASctl::omp_set_num_threads(nthreads)
  }
  cor.pearson <- function(x)
  {
    x = x - rowMeans(x)
    x = x / sqrt(rowSums(x^2))
    x.cor = tcrossprod(x)
    return(x.cor)
  }
  x <- as.matrix(x)
  if(!is.matrix(x)){
    warning("x is not like a matrix")
    return(NULL)
  }
  if(method=="pearson"){
    return(cor.pearson(x))
  }else if(method=="spearman"){
    return(cor.pearson(t(apply(x, 1, data.table::frank, na.last="keep"))))
  }else{
    warning("method must be pearson or spearman")
    return(NULL)
  }
}

#' dispaly message with time stamp
#' @param msg characters; message to display
loginfo <- function(msg) {
  timestamp <- sprintf("%s", Sys.time())
  msg <- paste0("[",timestamp, "] ", msg,"\n")
  cat(msg)
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


####### differential expression
#' differential expression analysis
#'
#' @importFrom plyr ldply
#' @importFrom stats aov TukeyHSD
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @param xdata data frame or matrix; rows for genes and columns for samples
#' @param xlabel factor; cluster label of the samples, with length equal to the number of columns in xdata
#' @param out.prefix character; if not NULL, write the result to the file(s). (default: NULL)
#' @param mod character;
#' @param F.FDR.THRESHOLD numeric; threshold of the adjusted p value of F-test. (default: 0.01)
#' @param HSD.FDR.THRESHOLD numeric; threshold of the adjusted p value of HSD-test (default: 0.01)
#' @param HSD.FC.THRESHOLD numeric; threshold of the absoute diff of HSD-test (default: 1)
#' @param verbose logical; whether output all genes' result. (default: F)
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
#' @param gid.mapping named character; gene id to gene symbol mapping. (default: NULL)
#' @return List with the following elements:
#' \item{aov.out}{data.frame, test result of all genes (rownames of xdata)}
#' \item{aov.out.sig}{format as aov.out, but only significant genes. }
findDEGenesByAOV <- function(xdata,xlabel,out.prefix=NULL,mod=NULL,
                             F.FDR.THRESHOLD=0.01,
                             HSD.FDR.THRESHOLD=0.01,
                             HSD.FC.THRESHOLD=1,
                             verbose=F,n.cores=NULL,
                             gid.mapping=NULL)
{
  clustNames <- unique(xlabel)
  xdata <- as.matrix(xdata)
  ### check rownames and colnames
  if(is.null(rownames(xdata))){
    stop("the xdata does not have rownames!!!")
  }
  if(length(table(xlabel))<2){
    return(NULL)
    #return(list(aov.out=ret.df,aov.out.sig=ret.df.sig))
  }
  if(!is.null(gid.mapping) && is.null(names(gid.mapping))) { names(gid.mapping)=gid.mapping }
  ###
  ## avoid conflict between threaded BLAs and foreach
  RhpcBLASctl::omp_set_num_threads(1)
  registerDoParallel(cores = n.cores)
  ret <- ldply(rownames(xdata),function(v){
    aov.out <- aov(y ~ g,data=data.frame(y=xdata[v,],g=xlabel))
    aov.out.s <- summary(aov.out)
    t.res.f <- unlist(aov.out.s[[1]]["g",c("F value","Pr(>F)")])
    aov.out.hsd <- TukeyHSD(aov.out)
    hsd.name <- rownames(aov.out.hsd$g)
    t.res.hsd <- c(aov.out.hsd$g[,"diff"],aov.out.hsd$g[,"p adj"])
    t.res.hsd.minP <- min(aov.out.hsd$g[,"p adj"])
    j <- which.min(aov.out.hsd$g[,"p adj"])
    t.res.hsd.minPDiff <- if(is.na(t.res.hsd.minP)) NaN else aov.out.hsd$g[j,"diff"]
    t.res.hsd.minPCmp <- if(is.na(t.res.hsd.minP)) NaN else rownames(aov.out.hsd$g)[j]
    ## whether cluster specific ?
    t.res.spe  <-  sapply(clustNames,function(v){ all( aov.out.hsd$g[grepl(v,hsd.name,perl=T),"p adj"] < HSD.FDR.THRESHOLD )  })
    ## wheter up across all comparison ?
    is.up <- sapply(clustNames,function(v){  all( aov.out.hsd$g[grepl(paste0(v,"-"),hsd.name),"diff"]>0 ) & all( aov.out.hsd$g[grepl(paste0("-",v),hsd.name),"diff"]<0 ) })
    is.down <- sapply(clustNames,function(v){  all( aov.out.hsd$g[grepl(paste0(v,"-"),hsd.name),"diff"]<0 ) & all( aov.out.hsd$g[grepl(paste0("-",v),hsd.name),"diff"]>0 ) })
    is.clusterSpecific <- (sum(t.res.spe,na.rm = T) == 1)
    if(is.clusterSpecific){
      t.res.spe.lable <- names(which(t.res.spe))
      if(is.up[t.res.spe.lable]) {
        t.res.spe.direction <- "UP"
      }else if(is.down[t.res.spe.lable]) {
        t.res.spe.direction <- "DOWN"
      }else{
        t.res.spe.direction <- "INCONSISTANT"
      }
    }else{
      t.res.spe.lable <- "NA"
      t.res.spe.direction <- "NA"
    }
    if(!is.null(mod) && mod=="cluster.specific") {
      structure(c(t.res.f,t.res.hsd,t.res.hsd.minP,t.res.hsd.minPDiff,t.res.hsd.minPCmp,
                  t.res.spe,is.clusterSpecific,t.res.spe.lable,t.res.spe.direction),
                names=c("F","F.pvalue",paste0("HSD.diff.",hsd.name),paste0("HSD.padj.",hsd.name),
                        "HSD.padj.min","HSD.padj.min.diff","HSD.padj.min.cmp",
                        paste0("cluster.specific.",clustNames),"is.clusterSpecific","cluster.lable","cluster.direction"))
    }else{
      structure(c(t.res.f,t.res.hsd,t.res.hsd.minP,t.res.hsd.minPDiff,t.res.hsd.minPCmp),
                names=c("F","F.pvalue",paste0("HSD.diff.",hsd.name),paste0("HSD.padj.",hsd.name),
                        "HSD.padj.min","HSD.padj.min.diff","HSD.padj.min.cmp"))
    }
  },.progress = "none",.parallel=T)
  #print(str(ret))

  if(!is.null(gid.mapping)){
    cnames <- gid.mapping[rownames(xdata)]
  }else{
    cnames <- rownames(xdata)
  }
  f.cnames.na <- which(is.na(cnames))
  cnames[f.cnames.na] <- rownames(xdata)[f.cnames.na]
  ret.df <- data.frame(geneID=rownames(xdata),geneSymbol=cnames,stringsAsFactors=F)
  ret.df <- cbind(ret.df,ret)
  rownames(ret.df) <- rownames(xdata)
  #print(str(ret.df))
  ## type conversion
  if(is.null(mod)){
    i <- 3:(ncol(ret.df)-1)
    ret.df[i]<-lapply(ret.df[i],as.numeric)
  }else{
    i <- 3:(ncol(ret.df)-(4+length(clustNames)))
    ret.df[i]<-lapply(ret.df[i],as.numeric)
    i <- (ncol(ret.df)-(length(clustNames)+2)):(ncol(ret.df)-2)
    ret.df[i]<-lapply(ret.df[i],as.logical)
  }
  ## adjust F test's p value
  ret.df$F.adjp <- 1
  ret.df.1 <- subset(ret.df,!is.na(F.pvalue))
  ret.df.1$F.adjp <- p.adjust(ret.df.1[,"F.pvalue"],method = "BH")
  ret.df.2 <- subset(ret.df,is.na(F.pvalue))
  ret.df <- rbind(ret.df.1,ret.df.2)
  ret.df <- ret.df[order(ret.df$F.adjp,-ret.df$F,ret.df$HSD.padj.min),]
  ### select
  ret.df.sig <- subset(ret.df,F.adjp<F.FDR.THRESHOLD & HSD.padj.min<HSD.FDR.THRESHOLD & abs(HSD.padj.min.diff)>=HSD.FC.THRESHOLD)
  ### output
  if(!is.null(out.prefix)){
    write.table(ret.df.sig,file = sprintf("%s.aov.sig.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    if(verbose){
      write.table(ret.df,file = sprintf("%s.aov.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    }
  }
  #print(str(ret.df))
  return(list(aov.out=ret.df,aov.out.sig=ret.df.sig))
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
    },error=function(e){ print("Error occur when using perplexity 5"); print(e); e })
  }
  return(ret)
}


#' Wraper for running SC3
#' @importFrom SC3 sc3 sc3_plot_consensus sc3_plot_silhouette sc3_plot_cluster_stability sc3_plot_markers
#' @importFrom plyr llply
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param out.prefix character, output prefix
#' @param n.cores integer, number of cors to use. (default: 8)
#' @param ks integer vector, number of clusters. (default: 2:10)
#' @param SC3.biology logical, SC3 parameter, whether calcualte biology. (default: T)
#' @param SC3.markerplot.width integer, SC3 parameter, with of the marker plot (default: 15)
#' @param verbose logical, whether verbose output. (default: F)
#' @details Run SC3 clustering pipeline
#' @return an object of \code{SingleCellExperiment} class with cluster labels and other info added.
#' @export
run.SC3 <- function(obj,assay.name="exprs",out.prefix=NULL,n.cores=8,ks=2:10,SC3.biology=T,SC3.markerplot.width=15,verbose=F)
{
  rownames.old <- rownames(obj)
  #### current SC3 need feature_symbol as rownames
  if(!"feature_symbol" %in% names(rowData(obj))){ rowData(obj)$feature_symbol <- rownames.old }
  rownames(obj) <- rowData(obj)$feature_symbol
  #### current SC3 use logcounts as dataset
  psu.logcounts <- F
  if(!"logcounts" %in% assayNames(obj)){
    assay(obj,"logcounts") <- assay(obj,assay.name)
    psu.logcounts <- T
  }
  #### current SC3 usde counts also
  psu.counts <- F
  if(!"counts" %in% assayNames(obj)){
    assay(obj,"counts") <- assay(obj,assay.name)
    psu.counts <- T
  }
  #### run
  obj <- sc3(obj, ks = ks, biology = SC3.biology, n_cores = n.cores,svm_max = 50000000,gene_filter = F)
  if(!is.null(out.prefix))
  {
    RhpcBLASctl::omp_set_num_threads(1)
    registerDoParallel(cores = n.cores)
    tryCatch({
      no.ret <- llply(ks,function(k){
        ###### sc3_plot_consensus is slow, not sure why
        png(sprintf("%s.consensus.k%d.png",out.prefix,k),width = 600,height = 480)
        sc3_plot_consensus(obj, k = k,  show_pdata = c( "sampleType", sprintf("sc3_%d_clusters",k),
                                                        sprintf("sc3_%s_log2_outlier_score",k)))
        dev.off()
        pdf(sprintf("%s.silhouette.k%d.pdf",out.prefix,k),width = 6,height = 6)
        sc3_plot_silhouette(obj, k = k)
        dev.off()
        p <- sc3_plot_cluster_stability(obj, k = k)
        ggsave(sprintf("%s.stability.k%d.pdf",out.prefix,k),width = 4,height = 3)
        if(SC3.biology){
          sc3_plot_markers(obj, k = k,auroc = 0.7,plot.extra.par = list(filename=sprintf("%s.markers.k%d.pdf",out.prefix,k),width=SC3.markerplot.width),
                           show_pdata = c( "sampleType",sprintf("sc3_%d_clusters",k), sprintf("sc3_%s_log2_outlier_score",k)))
        }
      },.progress = "none",.parallel=T)
    },error=function(e){
      cat(sprintf("Error occur in llply(ks,...).\n"))
      print(e)
    })
    if(verbose){ save(obj,file=sprintf("%s.verbose.sce.RData",out.prefix)) }
  }
  rownames(obj) <- rownames.old
  #### current SC3 use logcounts as dataset
  if(psu.logcounts){ assay(obj,"logcounts") <- NULL }
  #### current SC3 usde counts also
  if(psu.counts){ assay(obj,"counts") <- NULL }
  return(obj)
}

#' Wraper for running ZinbWave
#' @importFrom zinbwave zinbFit
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom BiocParallel MulticoreParam
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay to use for select genes (default: "exprs")
#' @param out.prefix character, output prefix
#' @param n.cores integer, number of cors to use. (default: 8)
#' @param zinbwave.K integer, zinbwave parameter, number of latent variables. (default: 20)
#' @param zinbwave.X character, zinbwave parameter, cell-level covariates. (default: "~patient")
#' @param verbose logical, whether verbose output. (default: F)
#' @details Run ZinbWave fitting
#' @return an object of class ZinbModel
#' @export
run.zinbWave <- function(obj,assay.name="exprs", vgene=NULL,out.prefix="./zinbwave",n.cores=8,
                         zinbwave.K=20,
                         zinbwave.X="~patient",verbose=F)
{
  if(is.null(vgene)){
    obj <- ssc.variableGene(obj,method = "sd",sd.n = 1500,assay.name = assay.name)
    vgene <- metadata(obj)$ssc$variable.gene$sd
  }
  RhpcBLASctl::omp_set_num_threads(1)
  #### fitting
  obj.zinb <- zinbFit(obj[vgene,], K=zinbwave.K, X=zinbwave.X, epsilon=1000,BPPARAM=MulticoreParam(n.cores),verbose=verbose)
  return(obj.zinb)
}




