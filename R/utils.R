#' @importFrom sscVis simple.removeBatchEffect
#' @export
sscVis::simple.removeBatchEffect

#' @importFrom sscVis run.cutree
#' @export
sscVis::run.cutree

#' @importFrom sscVis run.cutreeDynamic
#' @export
sscVis::run.cutreeDynamic

#' @importFrom sscVis loginfo
#' @export
sscVis::loginfo

#' @importFrom sscVis cor.BLAS
#' @export
sscVis::cor.BLAS



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
#' @importFrom utils write.table
#' @importFrom sscVis loginfo
#' @importFrom Matrix rowSums
#' @param xdata data frame or matrix; rows for genes and columns for samples
#' @param xlabel factor; cluster label of the samples, with length equal to the number of columns in xdata
#' @param batch factor; covariate. (default: NULL)
#' @param out.prefix character; if not NULL, write the result to the file(s). (default: NULL)
#' @param pmod character;
#' @param F.FDR.THRESHOLD numeric; threshold of the adjusted p value of F-test. (default: 0.01)
#' @param HSD.FDR.THRESHOLD numeric; threshold of the adjusted p value of HSD-test (default: 0.01)
#' @param HSD.FC.THRESHOLD numeric; threshold of the absoute diff of HSD-test (default: 1)
#' @param use.Kruskal logical; whether use Kruskal test for ranking genes (default: FALSE)
#' @param F.only logical; only perform F-test (default: FALSE)
#' @param verbose logical; whether output all genes' result. (default: F)
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param gid.mapping named character; gene id to gene symbol mapping. (default: NULL)
#' @return List with the following elements:
#' \item{aov.out}{data.frame, test result of all genes (rownames of xdata)}
#' \item{aov.out.sig}{format as aov.out, but only significant genes. }
findDEGenesByAOV <- function(xdata,xlabel,batch=NULL,out.prefix=NULL,pmod=NULL,
                             F.FDR.THRESHOLD=0.01,
                             HSD.FDR.THRESHOLD=0.01,
                             HSD.FC.THRESHOLD=1,
                             use.Kruskal=F,
                             F.only=F,
                             verbose=F,n.cores=NULL,
                             ncell.downsample=NULL,
                             gid.mapping=NULL)
{

  set.seed(9999)
  if(is.null(names(xlabel))) { names(xlabel) <- colnames(xdata) }
  if(!is.null(batch) && is.null(names(batch))) { names(batch) <- colnames(batch) }

  loginfo(sprintf("dim of xdata (input for findDEGenesByAOV): %d x %d",nrow(xdata),ncol(xdata)))

  #### downsample
  if(!is.null(ncell.downsample)){
    grp.list <- unique(as.character(xlabel))
    f.cell <- unname(unlist(sapply(grp.list, function(x) {
                x <- names(xlabel[xlabel == x])
                sample(x, min(length(x), ncell.downsample))
            })))
    xdata <- xdata[,f.cell,drop=F]
    xlabel <- xlabel[f.cell]
    batch <- batch[f.cell]
    loginfo(sprintf("dim of xdata (downsampled): %d x %d",nrow(xdata),ncol(xdata)))
  }  
  #### if xdata is dgCMatrix (from Matrix), Matrix::rowSums is required
  f.rm <- rowSums(xdata==0)==ncol(xdata)
  ####f.rm <- rowSds(xdata)==0
  xdata <- xdata[!f.rm,,drop=F]
  loginfo(sprintf("dim of xdata (remove all-zeros genes): %d x %d",nrow(xdata),ncol(xdata)))

  ####  
  clustNames <- unique(as.character(xlabel))
  ###### if xdata is data.frame, xdata[v,] will not be a vector
  if(is.data.frame(xdata)){ xdata <- as.matrix(xdata) }
  ### check rownames and colnames
  if(is.null(rownames(xdata))){
    stop("the xdata does not have rownames!!!")
  }
  if(length(table(xlabel))<2){
    return(NULL)
  }
  if(!is.null(gid.mapping) && is.null(names(gid.mapping))) { names(gid.mapping)=gid.mapping }
  ###
  ## avoid conflict between threaded BLAs and foreach
  RhpcBLASctl::omp_set_num_threads(1)
  registerDoParallel(cores = n.cores)
  ret <- as.data.table(ldply(rownames(xdata),function(v){

    xdata.v <- data.table(y=xdata[v,],g=xlabel,b=batch)
    xdata.v[,y.rnk:=rank(y)]

    t.res.kruskal <- kruskal.test(y ~ g,data=xdata.v)
    t.res.kruskal.p <- t.res.kruskal$p.value

    if(is.null(batch)){
        aov.out <- aov(y ~ g,data=xdata.v)
    }else{
        aov.out <- aov(y ~ g+b,data=xdata.v)
    }
    aov.out.s <- summary(aov.out)
    t.res.f <- unlist(aov.out.s[[1]]["g",c("F value","Pr(>F)")])

    dat.ret <- NULL
    dat.ret <- data.table("geneID"=v,"F"=t.res.f[1],"F.pvalue"=t.res.f[2],"Kruskal.pvalue"=t.res.kruskal.p)
    if(F.only){ return(dat.ret) }

    pairwise.cmp.tb <- NULL
    posthoc.pc.name <- NULL
    if(use.Kruskal){
        ###t.res.kruskalmc <- pgirmess::kruskalmc(y~g, data=data.frame(y=xdata[v,],g=xlabel),
        ###				       probs = 0.05, cont=NULL)
        #t.res.kruskal <- with(xdata.v,agricolae::kruskal(y,g,alpha=0.01,p.adj="BH",group=F))
        #t.res.kruskal.p <- t.res.kruskal$statistics$p.chisq

        #t.res.kruskal <- kruskal.test(y ~ g,data=xdata.v)
        #t.res.kruskal.p <- t.res.kruskal$p.value

        t.res.posthoc.raw <- with(xdata.v,agricolae::kruskal(y,g,alpha=HSD.FDR.THRESHOLD,p.adj="BH",group=F))
        posthoc.pc.name <- gsub(" ","",rownames(t.res.posthoc.raw$comparison))
        cmp.tb <- as.data.table(t.res.posthoc.raw$comparison)
        cmp.tb$cmp <- posthoc.pc.name
        cmp.tb <- cmp.tb[order(pvalue,-abs(Difference)),]
        pairwise.cmp.tb <- cmp.tb
        pairwise.cmp.tb$diff <- pairwise.cmp.tb$Difference
        pairwise.cmp.tb$p.adj <- pairwise.cmp.tb$pvalue
        ### pvalue actually is adjusted p value
        t.res.posthoc <- c(t.res.posthoc.raw$comparison[,"Difference"],t.res.posthoc.raw$comparison[,"pvalue"])
        t.res.posthoc.minP <- cmp.tb$pvalue[1]
        t.res.posthoc.minPDiff <- cmp.tb$Difference[1]
        t.res.posthoc.minPCmp <- cmp.tb$cmp[1]
        xdata.grp.mean <- xdata.v[,.(rnk.mean=mean(y.rnk)),by="g"][order(-rnk.mean),]
        cmp.minDiff <- cmp.tb[grepl(xdata.grp.mean$g[1],cmp,perl=T) &
                      grepl(xdata.grp.mean$g[2],cmp,perl=T) ,]
        t.res.posthoc.minDiff.P <- cmp.minDiff$pvalue[1]
        t.res.posthoc.minDiff <- cmp.minDiff$Difference[1]
        t.res.posthoc.minDiff.cmp <- cmp.minDiff$cmp[1]
    }else{
        aov.out.hsd <- TukeyHSD(aov.out)
        posthoc.pc.name <- rownames(aov.out.hsd$g)
        cmp.tb <- as.data.table(aov.out.hsd$g)
        cmp.tb$cmp <- posthoc.pc.name
        colnames(cmp.tb) <- make.names(colnames(cmp.tb))
        cmp.tb <- cmp.tb[order(p.adj,-abs(diff)),]
        pairwise.cmp.tb <- cmp.tb
        t.res.posthoc <- c(aov.out.hsd$g[,"diff"],aov.out.hsd$g[,"p adj"])
        t.res.posthoc.minP <- cmp.tb$p.adj[1]
        t.res.posthoc.minPDiff <- cmp.tb$diff[1]
        t.res.posthoc.minPCmp <- cmp.tb$cmp[1]
        xdata.grp.mean <- xdata.v[,.(mean=mean(y)),by="g"][order(-mean),]
        cmp.minDiff <- cmp.tb[grepl(xdata.grp.mean$g[1],cmp,perl=T) &
                      grepl(xdata.grp.mean$g[2],cmp,perl=T) ,]
        t.res.posthoc.minDiff.P <- cmp.minDiff$p.adj[1]
        t.res.posthoc.minDiff <- cmp.minDiff$diff[1]
        t.res.posthoc.minDiff.cmp <- cmp.minDiff$cmp[1]
    }
    {
        ## whether cluster specific ?
        t.res.spe  <-  sapply(clustNames,function(v){
                      all(pairwise.cmp.tb[grepl(v,cmp,perl=T),][["p.adj"]] < HSD.FDR.THRESHOLD) })
                      ##all( aov.out.hsd$g[grepl(v,hsd.name,perl=T),"p adj"] < HSD.FDR.THRESHOLD )
        ## wheter up across all comparison ?
        is.up <- sapply(clustNames,function(v){
                    all(pairwise.cmp.tb[grepl(paste0(v,"-"),posthoc.pc.name),][["diff"]] > 0) &
                    all(pairwise.cmp.tb[grepl(paste0("-",v),posthoc.pc.name),][["diff"]] < 0) })

        is.down <- sapply(clustNames,function(v){
                    all(pairwise.cmp.tb[grepl(paste0(v,"-"),posthoc.pc.name),][["diff"]] < 0) &
                    all(pairwise.cmp.tb[grepl(paste0("-",v),posthoc.pc.name),][["diff"]] > 0) })
			    
    }

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

    #dat.ret <- NULL
    #dat.ret <- data.table("geneID"=v,"F"=t.res.f[1],"F.pvalue"=t.res.f[2],"Kruskal.pvalue"=t.res.kruskal.p)
    dat.ret <- cbind(dat.ret, t(data.frame(structure(t.res.posthoc,names=c(paste0("PostHoc.diff.",posthoc.pc.name),
							     paste0("PostHoc.padj.",posthoc.pc.name))))))
    dat.ret <- cbind(dat.ret,data.table("PostHoc.padj.min"=t.res.posthoc.minP,
					"PostHoc.padj.min.diff"=t.res.posthoc.minPDiff,
					"PostHoc.padj.min.cmp"=t.res.posthoc.minPCmp,
					"PostHoc.minDiff.padj"=t.res.posthoc.minDiff.P,
					"PostHoc.minDiff"=t.res.posthoc.minDiff,
					"PostHoc.minDiff.cmp"=t.res.posthoc.minDiff.cmp))
    if(!is.null(pmod) && pmod=="cluster.specific") {
        dat.ret <- cbind(dat.ret,t(data.frame(structure(t.res.spe,names=paste0("cluster.specific.",clustNames)))))
        dat.ret <- cbind(dat.ret,data.table("is.clusterSpecific"=is.clusterSpecific,
                            "cluster.lable"=t.res.spe.lable,
                            "cluster.direction"=t.res.spe.direction))
    }
    
    return(dat.ret)
  },.progress = "none",.parallel=T))
  #print(str(ret))

  if(!is.null(gid.mapping)){
    cnames <- gid.mapping[rownames(xdata)]
  }else{
    cnames <- rownames(xdata)
  }
  f.cnames.na <- which(is.na(cnames))
  cnames[f.cnames.na] <- rownames(xdata)[f.cnames.na]
  
  ret.df <- data.table(geneID=rownames(xdata),geneSymbol=cnames)
  ret.df <- merge(ret.df,ret,by="geneID")
  
  ## adjust F test's p value
  ret.df$F.adjp <- 1
  ret.df$Kruskal.adjp <- 1
  ret.df.1 <- ret.df[!is.na(F.pvalue),]
  ret.df.1$F.adjp <- p.adjust(ret.df.1$F.pvalue,method = "BH")
  ret.df.2 <- ret.df[is.na(F.pvalue),]
  ret.df <- rbind(ret.df.1,ret.df.2)
  ##
  ret.df.3 <- ret.df[!is.na(Kruskal.pvalue),]
  ret.df.3$Kruskal.adjp <- p.adjust(ret.df.3$Kruskal.pvalue,method="BH")
  ret.df.4 <- ret.df[is.na(Kruskal.pvalue),]
  ret.df <- rbind(ret.df.3,ret.df.4)
  ##
  if(!F.only){
      ret.df <- ret.df[order(F.adjp,-F,PostHoc.padj.min),]
      ### select
      if(use.Kruskal){
          ret.df <- ret.df[order(Kruskal.adjp,-F,PostHoc.padj.min),]
          ret.df.sig <- ret.df[Kruskal.adjp<F.FDR.THRESHOLD &
                   PostHoc.padj.min<HSD.FDR.THRESHOLD &
                   abs(PostHoc.padj.min.diff)>=HSD.FC.THRESHOLD,]
      }else{
          ret.df.sig <- ret.df[F.adjp<F.FDR.THRESHOLD &
                   PostHoc.padj.min<HSD.FDR.THRESHOLD &
                   abs(PostHoc.padj.min.diff)>=HSD.FC.THRESHOLD,]
      }
  }else{
      ret.df <- ret.df[order(F.adjp,-F),]
      ret.df.sig <- ret.df[F.adjp < F.FDR.THRESHOLD,]
  }
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

#' calculate the AUC of one gene, using it as a classifier. code from SC3
#' @importFrom ROCR prediction performance
#' @importFrom stats aggregate wilcox.test
#' @param gene numeric; expression profile of one gene across samples
#' @param labels character; clusters of the samples belong to
#' @param use.rank logical; using the expression value itself or convert to rank value. (default: TRUE)
#' @param geneID character; geneID
#' @param method.Max character; method to find highest group, one of "mean", "median", "rank.mean" (default: mean)
getAUC <- function(gene, labels,use.rank=T,geneID="GeneX",method.Max="mean")
{
    requireNamespace("ROCR")

    if(use.rank){
        score <- rank(gene)
    }else{
        score <- gene
    }
    # Get average score for each cluster
    if(method.Max=="rank.mean"){
	score.rnk <- rank(gene)
	ms <- aggregate(score.rnk ~ labels, FUN = mean)
    }else{
	ms <- aggregate(score ~ labels, FUN = if(method.Max=="median") median else mean)
    }
    # Get cluster with highest average score
    posgroup <- ms[ms$score == max(ms$score), ]$labels
    # Return negatives if there is a tie for cluster with highest average score
    # (by definition this is not cluster specific)
    if(length(posgroup) > 1) {
        return (data.frame("geneID"=geneID,"AUC"=-1,"cluster"=posgroup[1],"score.p.value"=1,stringsAsFactors=F))
    }
    # Create 1/0 vector of truths for predictions, cluster with highest
    # average score vs everything else
    truth <- as.numeric(labels == posgroup)
    #Make predictions & get auc using RCOR package.
    pred <- ROCR::prediction(score,truth)
    val <- unlist(ROCR::performance(pred,"auc")@y.values)
    pval <- suppressWarnings(wilcox.test(score[truth == 1],
                                         score[truth == 0])$p.value)
    return(data.frame("geneID"=geneID,"AUC"=val,"cluster"=posgroup,"score.p.value"=pval,stringsAsFactors=F))
    #return(c(val,posgroup,pval))
}

#' For each gene, calculate the frequency of cells in each clusters are expressed.
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply
#' @param exp.bin numeric; binarized expression matrix, rows for genes and columns for samples. value 1 means expressed.
#' @param group character; clusters of the samples belong to
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
expressedFraction <- function(exp.bin,group,n.cores=NULL){
    requireNamespace("plyr")
    requireNamespace("doParallel")

    RhpcBLASctl::omp_set_num_threads(1)
    registerDoParallel(cores = n.cores)
    out.res <- ldply(rownames(exp.bin),function(v){
                .res <- aggregate(exp.bin[v,],by=list(group),FUN=function(x){ sum(x==1)/length(x) })
                structure(.res[,2],names=.res[,1])
			},.progress = "none",.parallel=T)
    rownames(out.res) <- rownames(exp.bin)
    colnames(out.res) <- sprintf("HiFrac.%s",colnames(out.res))
    out.res <- as.matrix(out.res)
    return(out.res)
}

#' For each gene, calculate the average expression of the expressor.
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr ldply
#' @param exp.bin numeric; binarized expression matrix, rows for genes and columns for samples. value 1 means expressed.
#' @param exp.norm numeric; expression matrix, rows for genes and columns for samples. original version of exp.bin.
#' @param group character; clusters of the samples belong to
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
expressedFraction.HiExpressorMean <- function(exp.bin,exp.norm,group,n.cores=NULL){
    requireNamespace("plyr")
    requireNamespace("doParallel")

    RhpcBLASctl::omp_set_num_threads(1)
    registerDoParallel(cores = n.cores)
    exp.bin[exp.bin<1] <- 0
    .exp <- exp.bin*exp.norm
    out.res <- ldply(rownames(.exp),function(v){
                .n <- aggregate(exp.bin[v,],by=list(group),FUN=function(x){ sum(x==1) })
                .res <- aggregate(.exp[v,],by=list(group),FUN=function(x){ sum(x) })
                structure(.res[,2]/.n[,2],names=.res[,1])
			},.progress = "none",.parallel=T)
    rownames(out.res) <- rownames(exp.bin)
    colnames(out.res) <- sprintf("AvgHiExpr.%s",colnames(out.res))
    out.res <- as.matrix(out.res)
    return(out.res)
}

#' Wraper for running random forest classifier (multiple core version)
#'
#' @importFrom randomForest randomForest
#' @importFrom parallel mclapply
#' @param ntree integer; Number of trees to grow.
#' @param ncore integer; Number of cores to use.
#' @param ... passed to randomForest::randomForest
mcrf <- function(ntree, ncore, ...) {
	ntrees <- rep(ntree%/%ncore, ncore) + c(rep(0, ncore-ntree%%ncore), rep(1, ntree%%ncore))
	rfs <- mclapply(ntrees, function(nn) {
#						a.mod <- randomForest::randomForest(ntree = nn, keep.inbag=T, ...)
#						pred <- predict(a.mod,x, predict.all=TRUE)
#						oob.data <- as.data.table(ldply(seq_len(nn),function(i){
#							tree_pred_i <- pred$individual[ , i]
#							oob_idx <- a.mod$inbag[ , i] == 0
#							data.frame(oob.pred=tree_pred_i[oob_idx],
#									   oob.actual=a.mod$y[oob_idx])
#						}))
#						conf_mtx <- oob.data[,table(oob.pred,oob.actual)]
#						error_rates <- 1 - diag(conf_mtx) / colSums(conf_mtx)
#						error_rate <- 1 - sum(diag(conf_mtx)) / sum(conf_mtx)
						a.mod <- randomForest::randomForest(ntree = nn, ...)
						return(a.mod)
	}, mc.cores = ncore)
	do.call(randomForest::combine, rfs)
}

####### classification functions

#' Wraper for running random forest classifier
#'
#' @importFrom varSelRF varSelRF
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @importFrom ranger ranger
#' @param xdata data frame or matrix; data used for training, with sample id in rows and variables in columns
#' @param xlabel factor; classification label of the samples, with length equal to the number of rows in xdata
#' @param ydata data frame or matrix; data to be predicted the label, same format as xdata
#' @param do.norm logical; whether perform Z score normalization on data
#' @param ntree integer; parameter of varSelRF::varSelRF
#' @param ntreeIterat integer; parameter of varSelRF::varSelRF
#' @param selectVar logical; wheter select variables (default: TRUE)
#' @param use.ranger logical; wheter use ranger (default: TRUE)
#' @param ncores intersect; numer of cores to use (default: NULL)
#' @param xdata.test data frame or matrix; data used for test, with sample id in rows and variables in columns
#' @param xlabel.test factor; classification label of the test samples, with length equal to the number of rows in xdata.test
#' @return List with the following elements:
#' \item{ylabel}{ppredicted labels of the samples in ydata}
#' \item{rfsel}{trained model; output of varSelRF()}
run.RF <- function(xdata, xlabel, ydata, do.norm=F,
				   ntree = 500, ntreeIterat = 200, selectVar=T,
				   ncores=NULL,use.ranger=T,
				   xdata.test=NULL,xlabel.test=NULL)
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
  pred.test <- NULL
  ### random forest
  if(selectVar){
	  rf.model <- varSelRF::varSelRF(xdata, xlabel,
								  ntree = ntree, ntreeIterat = ntreeIterat,
								  whole.range = FALSE,keep.forest = T)
	  #rfsel$selected.vars %>% str %>% print
	  #rfsel$initialImportances %>% head %>% print
	  #rfsel$rf.model$confusion %>% print
	  yres <- predict(rf.model$rf.model, newdata = ydata[,rf.model$selected.vars],type = "prob")
  }else{
	  if(!is.null(ncores) && ncores > 1){
		  if(use.ranger){
			  rf.model <- ranger::ranger(x=xdata,y=xlabel,probability=T,num.trees=ntree,
										 importance="impurity",num.threads=ncores,classification=T)
			  yres <- predict(rf.model,ydata,num.threads=ncores)$predictions
			  rf.model$importance <- data.frame(x=rf.model$variable.importance)
			  colnames(rf.model$importance) <- rf.model$importance.mode
			  if(!is.null(xdata.test)){
				  pred.test <- predict(rf.model,xdata.test,num.threads=ncores)$predictions
			  }
		  }else{
			  rf.model <- mcrf(ntree, ncores, x=xdata,y=xlabel,importance=T,proximity=T)
			  yres <- predict(rf.model,ydata,type="prob")
			  if(!is.null(xdata.test)){
				  pred.test <- predict(rf.model,xdata.test,type="prob")
			  }
		  }
	  }else{
		  rf.model <- randomForest::randomForest(xdata, xlabel,importance=T,proximity=T,ntree=ntree)
		  yres <- predict(rf.model,ydata,type="prob")
	  }
  }
  cls.set <- colnames(yres)
  ylabel <- apply(yres,1,function(x){ cls.set[which.max(x)] })
  names(ylabel) <- rownames(ydata)
  ret.list <- list("ylabel"=ylabel,"rfsel"=rf.model,"yres"=yres)
  if(!is.null(pred.test)){
	  pred.test.res <- data.table(pred=apply(pred.test,1,function(x){ as.character(cls.set[which.max(x)]) }),
								  actual=as.character(xlabel.test),
								  group=as.character(xlabel.test))
	  conf.mtx.test <- pred.test.res[,table(pred,actual)]
	  error.rate.perCls <- pred.test.res[,.(N=.N, error.rate=sum(.SD$pred!=.SD$actual)/.N),by="group"]
	  error.rate.overall <- sum(pred.test.res[,pred!=actual])/nrow(pred.test.res)
	  ret.list[["conf.mtx"]] <- conf.mtx.test
	  ret.list[["error.rate.perCls"]] <- error.rate.perCls
	  ret.list[["error.rate.overall"]] <- error.rate.overall
  }
  return(ret.list)
}

#' Wraper for running svm
#'
#' @importFrom e1071 svm
#' @importFrom stats predict
#' @param xdata data frame or matrix; data used for training, with sample id in rows and variables in columns
#' @param xlabel factor; classification label of the samples, with length equal to the number of rows in xdata
#' @param ydata data frame or matrix; data to be predicted the label, same format as xdata
#' @param kern character; which kernel to use, can be one of linear, polynomial, radial and sigmoid (default: "linear")
#' @return List with the following elements:
#' \item{ylabel}{ppredicted labels of the samples in ydata}
#' \item{rfsel}{trained model; output of varSelRF()}
run.SVM <- function(xdata, xlabel, ydata,kern="linear")
{
  f.g <- intersect(colnames(xdata),colnames(ydata))
  xdata <- xdata[,f.g,drop=F]
  ydata <- ydata[,f.g,drop=F]
  model <- e1071::svm(xdata, xlabel, kernel = kern)
  ylabel <- predict(model, newdata=ydata)
  names(ylabel) <- rownames(ydata)
  return(list("ylabel"=ylabel,"svm"=model))
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
#' @importFrom stats prcomp
#' @param idata matrix; expression data with sample id in rows and variables in columns
#' @param tSNE.usePCA whether perform PCA before tSNE (default: T)
#' @param tSNE.perplexity perplexity parameter of tSNE (default: 30)
#' @param method method to be used. one of "Rtsne" and "FIt-SNE" (default: "Rtsne")
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
#' @param out.prefix character; output prefix (default: NULL)
#' @param ... parameters passed to Rtsne or fftRtsne
#' @return If successful same as the return value of Rtsne(); otherwise NULL
run.tSNE <- function(idata,tSNE.usePCA=T,tSNE.perplexity=30,method="Rtsne",n.cores=NULL,out.prefix=NULL,...){
  ret <- NULL
  if(is.null(n.cores)){ n.cores <- 1 }
  if(method=="Rtsne"){
      tryCatch({
        ret <- Rtsne::Rtsne(idata, pca = tSNE.usePCA, num_threads=n.cores, 
							perplexity = tSNE.perplexity,...)$Y
      },error=function(e){
        #cat("Perplexity is too large; try to use smaller perplexity 5\n")
      })
      if(is.null(ret)){
		cat(sprintf("Warning: re-run Rtsne() using perplexity=5\n"))
        tryCatch({
          ret <- Rtsne::Rtsne(idata, pca = tSNE.usePCA, num_threads=n.cores, perplexity = 5,...)$Y
        },error=function(e){ print("Error occur when using perplexity 5"); print(e); e })
      }
  }else if(method=="FIt-SNE"){
      if(tSNE.usePCA){
        pca.res <- prcomp(idata)
        pca.npc <- min(50,ncol(pca.res$x))
        X <- pca.res$x[,1:pca.npc,drop=F]
      }else{
        X <- idata
      }
      tryCatch({
		ret <- fftRtsne(X,perplexity=tSNE.perplexity,nthreads=n.cores,out_prefix=out.prefix,...)
	  },error=function(e){
	  })
	  if(is.null(ret)){
		cat(sprintf("Warning: re-run fftRtsne() using perplexity=5\n"))
        tryCatch({
		  ret <- fftRtsne(X,perplexity=5,nthreads=n.cores,out_prefix=out.prefix,...)
        },error=function(e){ print("Error occur when using perplexity 5"); print(e); e })
      }
  }
  return(ret)
}


#' Wraper for running SC3
#' @importFrom SC3 sc3 sc3_plot_consensus sc3_plot_silhouette sc3_plot_cluster_stability sc3_plot_markers
#' @importFrom plyr llply
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom SummarizedExperiment `rowData<-` `assay<-`
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
#          sc3_plot_markers(obj, k = k,auroc = 0.7,plot.extra.par = list(filename=sprintf("%s.markers.k%d.pdf",out.prefix,k),
#                                                                        width=SC3.markerplot.width),
#                           show_pdata = c( "sampleType",sprintf("sc3_%d_clusters",k), sprintf("sc3_%s_log2_outlier_score",k)))
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


#' Wraper for running FIt-SNE. Code from KlugerLab (https://github.com/KlugerLab/FIt-SNE)
#' @param X matrix; samples in rows and variables in columns
#' @param dims integer; dimentionality of the returned matrix
#' @param perplexity double; perplexity parameter of tSNE (effective nearest neighbours)
#' @param theta double; theta
#' @param max_iter integer; max_iter
#' @param fft_not_bh logical; fft_not_bh (default: TRUE)
#' @param ann_not_vptree logical; (default: TRUE)
#' @param stop_early_exag_iter integer; (default 250)
#' @param exaggeration_factor double;
#' @param no_momentum_during_exag logical;
#' @param start_late_exag_iter double;
#' @param late_exag_coeff double;
#' @param mom_switch_iter double;
#' @param momentum double;
#' @param final_momentum double;
#' @param learning_rate double;
#' @param n_trees integer;
#' @param search_k double;
#' @param rand_seed double;
#' @param nterms integer;
#' @param intervals_per_integer integer;
#' @param min_num_intervals integer;
#' @param K integer;
#' @param sigma double;
#' @param initialization matrix
#' @param out_prefix character; temporary files prefix for fast_tsne (default: NULL)
#' @param load_affinities character;
#' @param fast_tsne_path character; full path of the installed fast_tsne programe (default: NULL)
#' @param nthreads integer; number of threads (default: 0)
#' @param perplexity_list list;
#' @param get_costs logical;
#' @param data_path character;
#' @param result_path character;
#' @param df numeric;
#' @importFrom utils file_test
#' @details Run FIt-SNE
#' @return a matrix with samples in rows and tSNE coordinate in columns
#' @export
fftRtsne <- function(X, 
		     dims = 2, perplexity = 30, theta = 0.5,
		     max_iter = 1000,
		     fft_not_bh = TRUE,
		     ann_not_vptree = TRUE,
		     stop_early_exag_iter = 250,
		     exaggeration_factor = 12.0, no_momentum_during_exag = FALSE,
		     start_late_exag_iter = -1.0, late_exag_coeff = 1.0,
         mom_switch_iter = 250, momentum = 0.5, final_momentum = 0.8, learning_rate = 200,
		     n_trees = 50, search_k = -1, rand_seed = -1,
		     nterms = 3, intervals_per_integer = 1, min_num_intervals = 50, 
		     K = -1, sigma = -30, initialization = NULL,
		     data_path = NULL, result_path = NULL,out_prefix = NULL,
		     load_affinities = NULL,
		     fast_tsne_path = NULL, nthreads = 0, perplexity_list = NULL, 
         get_costs = FALSE, df = 1.0) {
  
  version_number <- '1.1.0'

	data_path <- tempfile(pattern=sprintf("%s.fftRtsne_data_",
										  if(is.null(out_prefix)) "" else out_prefix),
                        fileext='.dat')
	result_path <- tempfile(pattern=sprintf("%s.fftRtsne_result_",
											if(is.null(out_prefix)) "" else out_prefix),
                          fileext='.dat')

	if (is.null(fast_tsne_path)) {
		fast_tsne_path <- system2('which', 'fast_tsne', stdout = TRUE)
	}
	fast_tsne_path <- normalizePath(fast_tsne_path)
	if (!file_test('-x', fast_tsne_path)) {
		stop(fast_tsne_path, " does not exist or is not executable; check your fast_tsne_path parameter")
	}

	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

	if (!is.numeric(theta) || (theta < 0.0) || (theta > 1.0) ) { stop("Incorrect theta.")}
	if (nrow(X) - 1 < 3 * perplexity) { stop("Perplexity is too large.")}
	if (!is.matrix(X)) { stop("Input X is not a matrix")}
	if (!(max_iter > 0)) { stop("Incorrect number of iterations.")}
	if (!is.wholenumber(stop_early_exag_iter) || stop_early_exag_iter < 0) { stop("stop_early_exag_iter should be a positive integer")}
	if (!is.numeric(exaggeration_factor)) { stop("exaggeration_factor should be numeric")}
	if (!is.numeric(df)) { stop("df should be numeric")}
	if (!is.wholenumber(dims) || dims <= 0) { stop("Incorrect dimensionality.")}
	if (search_k == -1) {
    if (perplexity > 0) {
      search_k <- n_trees * perplexity * 3
    } else if (perplexity == 0) {
      search_k <- n_trees * max(perplexity_list) * 3
    } else { 
      search_k <- n_trees * K
    }
  }

	if (fft_not_bh) {
	  nbody_algo <- 2
	} else {
	  nbody_algo <- 1
	}

	if (is.null(load_affinities)) {
		load_affinities <- 0
	} else {
		if (load_affinities == 'load') {
			load_affinities <- 1
		} else if (load_affinities == 'save') {
			load_affinities <- 2
		} else {
			load_affinities <- 0
		}
	}
	
	if (ann_not_vptree) {
	  knn_algo <- 1
	} else {
	  knn_algo <- 2
	}
	tX <- as.numeric(t(X))

	f <- file(data_path, "wb")
	n <- nrow(X)
	D <- ncol(X)
	writeBin(as.integer(n), f, size = 4)
	writeBin(as.integer(D), f, size = 4)
	writeBin(as.numeric(theta), f, size = 8) #theta
	writeBin(as.numeric(perplexity), f, size = 8)

  if (perplexity == 0) {
  	writeBin(as.integer(length(perplexity_list)), f, size = 4)
    writeBin(perplexity_list, f) 
  }

	writeBin(as.integer(dims), f, size = 4)
	writeBin(as.integer(max_iter), f, size = 4)
	writeBin(as.integer(stop_early_exag_iter), f, size = 4)
	writeBin(as.integer(mom_switch_iter), f, size = 4)
	writeBin(as.numeric(momentum), f, size = 8)
	writeBin(as.numeric(final_momentum), f, size = 8)
	writeBin(as.numeric(learning_rate), f, size = 8)
	writeBin(as.integer(K), f, size = 4) #K
	writeBin(as.numeric(sigma), f, size = 8) #sigma
	writeBin(as.integer(nbody_algo), f, size = 4)  #not barnes hut
	writeBin(as.integer(knn_algo), f, size = 4) 
	writeBin(as.numeric(exaggeration_factor), f, size = 8) #compexag
	writeBin(as.integer(no_momentum_during_exag), f, size = 4) 
	writeBin(as.integer(n_trees), f, size = 4) 
	writeBin(as.integer(search_k), f, size = 4) 
	writeBin(as.integer(start_late_exag_iter), f, size = 4) 
	writeBin(as.numeric(late_exag_coeff), f, size = 8) 
	
	writeBin(as.integer(nterms), f, size = 4) 
	writeBin(as.numeric(intervals_per_integer), f, size = 8) 
	writeBin(as.integer(min_num_intervals), f, size = 4) 
	writeBin(tX, f) 
	writeBin(as.integer(rand_seed), f, size = 4) 
  writeBin(as.numeric(df), f, size = 8)
	writeBin(as.integer(load_affinities), f, size = 4) 
	if (!is.null(initialization)) { writeBin( c(t(initialization)), f) }		
	close(f) 

	flag <- system2(command = fast_tsne_path, 
	                args = c(version_number, data_path, result_path, nthreads))
	if (flag != 0) {
		stop('tsne call failed')
	}
	f <- file(result_path, "rb")
	n <- readBin(f, integer(), n = 1, size = 4)
	d <- readBin(f, integer(), n = 1, size = 4)
	Y <- readBin(f, numeric(), n = n * d)
  Y <- t(matrix(Y, nrow = d))
  if (get_costs) {
    readBin(f, integer(), n = 1, size = 4)
    costs <- readBin(f, numeric(), n = max_iter, size = 8)
    Yout <- list(Y = Y, costs = costs)
  } else {
    Yout <- Y
  }
  close(f)
  file.remove(data_path)
  file.remove(result_path)
  Yout
}

#' @importFrom sscVis simple.removeBatchEffect
#' @export
sscVis::simple.removeBatchEffect



#' run limma, given an expression matrix
#' @param xdata data frame or matrix; rows for genes and columns for samples
#' @param xlabel factor; cluster label of the samples, with length equal to the number of columns in xdata
#' @param batch factor; covariate. (default: NULL)
#' @param out.prefix character; if not NULL, write the result to the file(s). (default: NULL)
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param T.fdr numeric; threshold of the adjusted p value of moderated t-test (default: 0.05)
#' @param T.logFC numeric; threshold of the absoute diff (default: 1)
#' @param T.expr numeric; threshold for binarizing exprs (default: 0.3)
#' @param T.bin.useZ logical; wheter use the z-score version of assay.namme for binarizing exprs (default: T)
#' @param verbose integer; verbose (default: 0)
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
#' @param group character; group of interest, if NULL the last group will be used (default: NULL)
#' @param gid.mapping named character; gene id to gene symbol mapping. (default: NULL)
#' @param do.voom logical; perform voom transfromation (default: FALSE)
#' @param rn.seed integer; random number seed (default: 9999)
#' @details diffeerentially expressed genes dectection using limma
#' @return a matrix with dimention as input ( samples in rows and variables in columns)
#' @importFrom data.table data.table as.data.table
#' @importFrom grDevices pdf dev.off
#' @importFrom limma lmFit eBayes topTable voom makeContrasts contrasts.fit
#' @importFrom stats model.matrix
#' @importFrom utils write.table
#' @importFrom matrixStats rowVars
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom BiocParallel MulticoreParam register multicoreWorkers
#' @export
run.limma.matrix <- function(xdata,xlabel,batch=NULL,out.prefix=NULL,ncell.downsample=NULL,
                             T.fdr=0.05,T.logFC=1,T.expr=0.3,T.bin.useZ=T,
							 verbose=0,n.cores=NULL,group=NULL,
                             gid.mapping=NULL, do.voom=F,rn.seed=9999)
{
	#suppressPackageStartupMessages(require("limma"))
	#suppressPackageStartupMessages(require("dplyr"))
	#suppressPackageStartupMessages(require("BiocParallel"))

    if(ncol(xdata)!=length(xlabel)){
        warning(sprintf("xdata and xlabel is not consistent!\n"))
        return(NULL)
    }
    names(xlabel) <- colnames(xdata)
    if(!is.null(ncell.downsample)){
        set.seed(rn.seed)
        grp.list <- unique(xlabel)
        f.cell <- unlist(lapply(grp.list,function(x){
                             x <- names(xlabel[xlabel==x])
                             sample(x,min(length(x),ncell.downsample)) }))
        xlabel <- xlabel[f.cell]
        if(!is.null(batch)){
            names(batch) <- colnames(xdata)
            batch <- batch[f.cell]
        }
        xdata <- xdata[,f.cell,drop=F]
    }
    xdata <- as.matrix(xdata)
    f.gene <- rowVars(xdata)>0 & apply(xdata,1,function(x){ !any(is.na(x)) })
    xdata <- xdata[f.gene,,drop=F]
    ####
	if("_control" %in% xlabel){
		x.levels <- c("_control",setdiff(unique(sort(xlabel)),"_control"))
	}else{
		x.levels <- unique(sort(xlabel))
	}
	group.dis <- NULL
    if(!is.null(group)){
		tmp.group <- unlist(strsplit(as.character(group),":"))
		group <- tmp.group[1]
		if(length(tmp.group)>1){
			group.dis <- tmp.group[2]
		}
        x.levels <- c(setdiff(x.levels,group),as.character(group))
    }

    if(is.null(batch)){
        design.df <- data.frame(cellID=colnames(xdata),
                                group=factor(xlabel,levels=x.levels,
                                             labels=c(sprintf("G%04d",seq_len(length(x.levels)-1)),"II")),
                                xlabel=xlabel,
                                stringsAsFactors=F)
	    design <- model.matrix(~group,data=design.df)
    }else{
        design.df <- data.frame(cellID=colnames(xdata),
                                group=factor(xlabel,levels=x.levels,
                                             labels=c(sprintf("G%04d",seq_len(length(x.levels)-1)),"II")),
                                batch=factor(batch,levels=unique(sort(batch)),
                                             labels=sprintf("B%04d",seq_along(unique(sort(batch))))),
                                xlabel=xlabel,
                                stringsAsFactors=F)
	    design <- model.matrix(~batch+group,data=design.df)
    }
    group.label <- x.levels[length(x.levels)]
	colnames(design)[ncol(design)]<-"II"
    colnames(design) <- gsub("[)(]","",colnames(design))

    RhpcBLASctl::omp_set_num_threads(1)
	register(MulticoreParam(if(!is.null(n.cores)) n.cores else multicoreWorkers()))
    if(do.voom){
        if(!is.null(out.prefix)){
            pdf(sprintf("%s.voom.pdf",out.prefix),width = 7,height = 7)
            v <- voom(xdata, design, plot=TRUE)
            dev.off()
        }
	    fit <- lmFit(v, design)
    }else{
	    fit <- lmFit(xdata, design)
    }
    fit.00 <- fit

    col.group <- grep("^group",colnames(design),value=T)
    if(length(x.levels)>2){
        contrast.str  <-  sprintf("II-(%s)/%d",
                                  paste(col.group,collapse="+"),
                                  length(col.group))
        contrasts.matrix <- makeContrasts(contrasts=contrast.str,levels=design)
        fit <- contrasts.fit(fit,contrasts=contrasts.matrix)
        fit <- eBayes(fit)
	    all.table  <- topTable(fit, number = Inf, sort.by = "p", p.value = 1)
    }else{
	    fit <- eBayes(fit)
	    all.table  <- topTable(fit, coef = "II", number = Inf, sort.by = "p", p.value = 1)
    }
	all.table <- cbind(data.table(geneID=rownames(all.table),
                                  geneSymbol=if(is.null(gid.mapping)) rownames(all.table) else gid.mapping[rownames(all.table)],
                                  cluster=if(is.null(group.dis)) group.label else group.dis,
                                  stringsAsFactors = F),
                       all.table)

    .getStatistics <- function(inDat,str.note=NULL,stat="mean"){
        .Grp.stat <- t(apply(inDat,1,function(x){
                                 .mexp <- aggregate(x ~ design.df$xlabel,
													FUN = if(stat=="mean") mean else if(stat=="sd") sd else if(stat=="length") length else if(stat=="freq") function(x){ sum(x>T.expr)/length(x) })
                                 structure(.mexp[,2],names=sprintf("%s.%s",stat,.mexp[,1]))
								  }))
        if(!is.null(str.note)){
            colnames(.Grp.stat) <- sprintf("%s.%s",colnames(.Grp.stat),str.note)
        }
        .Grp.stat.df <- data.table(geneID=rownames(.Grp.stat),stringsAsFactors = F)
        .Grp.stat.df <- cbind(.Grp.stat.df,.Grp.stat)
    }
    if(verbose>0){
        if(!is.null(batch)){
            betaV <- fit.00$coefficients
            betaV[is.na(betaV)] <- 0
            idx.batch <- grep("^batch",colnames(design),value = T)
            xdata.rmBE <- as.matrix(xdata) - betaV[,idx.batch,drop=F] %*% t(design[,idx.batch,drop=F])
            .Grp.mean.df <- .getStatistics(xdata.rmBE,str.note="rmBE")
			xdata.scale <- t(scale(t(xdata.rmBE)))
            .Grp.mean.scale.df <- .getStatistics(xdata.scale,str.note="scale")
            all.table <- merge(all.table,.Grp.mean.df)
            all.table <- merge(all.table,.Grp.mean.scale.df)
			all.table$meanExp <- all.table[[sprintf("mean.%s.rmBE",group.label)]]
			all.table$meanScale <- all.table[[sprintf("mean.%s.scale",group.label)]]
			if(verbose>1){
				.Grp.sd.df <- .getStatistics(xdata.rmBE,str.note="rmBE",stat="sd")
				all.table <- merge(all.table,.Grp.sd.df)
				colName.sd <- grep("^sd.",colnames(all.table),value=T)
				all.table$SNR <- all.table$logFC/rowSums(all.table[,colName.sd,with=F])
				.Grp.n.df <- .getStatistics(xdata.rmBE,stat="length")
				all.table <- merge(all.table,.Grp.n.df)
				if(T.bin.useZ){
					.Grp.freq.df <- .getStatistics(xdata.scale,stat="freq")
				}else{
					.Grp.freq.df <- .getStatistics(xdata,stat="freq")
				}
				all.table <- merge(all.table,.Grp.freq.df)
			}
        }else{
            .Grp.mean.df <- .getStatistics(xdata)
			xdata.scale <- t(scale(t(xdata)))
            .Grp.mean.scale.df <- .getStatistics(xdata.scale,str.note="scale")
            all.table <- merge(all.table,.Grp.mean.df)
            all.table <- merge(all.table,.Grp.mean.scale.df)
			all.table$meanExp <- all.table[[sprintf("mean.%s",group.label)]]
			all.table$meanScale <- all.table[[sprintf("mean.%s.scale",group.label)]]
			if(verbose>1){
				.Grp.sd.df <- .getStatistics(xdata,stat="sd")
				all.table <- merge(all.table,.Grp.sd.df)
				colName.sd <- grep("^sd.",colnames(all.table),value=T)
				all.table$SNR <- all.table$logFC/rowSums(all.table[,colName.sd,with=F])
				.Grp.n.df <- .getStatistics(xdata,stat="length")
				all.table <- merge(all.table,.Grp.n.df)
				if(T.bin.useZ==T){
					.Grp.freq.df <- .getStatistics(xdata.scale,stat="freq")
				}else{
					.Grp.freq.df <- .getStatistics(xdata,stat="freq")
				}
				all.table <- merge(all.table,.Grp.freq.df)
			}
        }
    }

    all.table <- all.table[order(adj.P.Val,-t,-logFC),]

    #print(head(all.table))
    sig.table <- all.table[adj.P.Val<T.fdr & abs(logFC)>T.logFC,]
    if(!is.null(out.prefix))
    {
		conn <- gzfile(sprintf("%s.limma.all.txt.gz",out.prefix),"w")
        write.table(all.table,conn,row.names = F,quote = F,sep = "\t")
		close(conn)
		conn <- gzfile(sprintf("%s.limma.sig.txt.gz",out.prefix),"w")
        write.table(sig.table,conn,row.names = F,quote = F,sep = "\t")
		close(conn)
#		if(verbose>1){
#			saveRDS(fit,file=sprintf("%s.fit.rds",out.prefix))
#		}
    }
	ret.dat <- list(all=all.table,sig=sig.table)
	if(verbose>1){ ret.dat$fit <- fit }
    return(ret.dat)
}

#' run DE given an expression matrix
#' @param xdata data frame or matrix; rows for genes and columns for samples
#' @param xlabel factor; cluster label of the samples, with length equal to the number of columns in xdata
#' @param batch factor; covariate. (default: NULL)
#' @param out.prefix character; if not NULL, write the result to the file(s). (default: NULL)
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param T.fdr numeric; threshold of the adjusted p value of moderated t-test (default: 0.05)
#' @param T.logFC numeric; threshold of the absoute diff (default: 1)
#' @param verbose integer; verbose (default: 0)
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
#' @param group character; group of interest, if NULL the last group will be used (default: NULL)
#' @param gid.mapping named character; gene id to gene symbol mapping. (default: NULL)
#' @param rn.seed integer; random number seed (default: 9999)
#' @param method character; . (default: "lm")
#' @details diffeerentially expressed genes dectection using naive methods
#' @return a matrix with dimention as input ( samples in rows and variables in columns)
#' @importFrom data.table data.table as.data.table
#' @importFrom grDevices pdf dev.off
#' @importFrom stats model.matrix
#' @importFrom utils write.table
#' @importFrom matrixStats rowVars
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom BiocParallel MulticoreParam register multicoreWorkers
#' @export
run.DE.matrix <- function(xdata,xlabel,batch=NULL,out.prefix=NULL,ncell.downsample=NULL,
			  T.fdr=0.05,T.logFC=1,T.expr=0.3,T.bin.useZ=T,
			  verbose=0,n.cores=NULL,group=NULL,
			  gid.mapping=NULL, do.voom=F,rn.seed=9999,method="lm")
{

    if(ncol(xdata)!=length(xlabel)){
	warning(sprintf("xdata and xlabel is not consistent!\n"))
	return(NULL)
    }
    names(xlabel) <- colnames(xdata)
    if(!is.null(ncell.downsample)){
	set.seed(rn.seed)
	grp.list <- unique(xlabel)
	f.cell <- unlist(lapply(grp.list,function(x){
			     x <- names(xlabel[xlabel==x])
			     sample(x,min(length(x),ncell.downsample)) }))
	xlabel <- xlabel[f.cell]
	if(!is.null(batch)){
	    names(batch) <- colnames(xdata)
	    batch <- batch[f.cell]
	}
	xdata <- xdata[,f.cell,drop=F]
    }
    xdata <- as.matrix(xdata)
    f.gene <- rowVars(xdata)>0 & apply(xdata,1,function(x){ !any(is.na(x)) })
    xdata <- xdata[f.gene,,drop=F]
    ####
    if("_control" %in% xlabel){
	x.levels <- c("_control",setdiff(unique(sort(xlabel)),"_control"))
    }else{
	x.levels <- unique(sort(xlabel))
    }
    ### group name for display purpose
    group.dis <- NULL
    if(!is.null(group)){
	tmp.group <- unlist(strsplit(as.character(group),":"))
	group <- tmp.group[1]
	if(length(tmp.group)>1){
	    group.dis <- tmp.group[2]
	}
	x.levels <- c(setdiff(x.levels,group),as.character(group))
    }

    group.label <- x.levels[length(x.levels)]

    RhpcBLASctl::omp_set_num_threads(1)
    register(MulticoreParam(if(!is.null(n.cores)) n.cores else multicoreWorkers()))

    {
	all.table <- as.data.table(ldply(rownames(xdata),function(v){
	    if(is.null(batch)){
		xdata.v <- data.table(y=xdata[v,],g=xlabel)
	    }else{
		xdata.v <- data.table(y=xdata[v,],g=xlabel,b=batch)
	    }
	    
	    {
		test.out <- wilcox.test(data=xdata.v,y~g)
		lm.out <- lm(data=xdata.v,y~g)
		lm.out.s <- summary(lm.out)
		data.table(geneID=v,
			   geneSymbol=if(is.null(gid.mapping)) v else gid.mapping[v],
			   cluster=if(is.null(group.dis)) group.label else group.dis,
			   logFC=lm.out.s$coefficients[sprintf("g%s",group),"Estimate"],
			   t=lm.out.s$coefficients[sprintf("g%s",group),"t value"],
			   p.value=lm.out.s$coefficients[sprintf("g%s",group),"Pr(>|t|)"],
			   wilcox.p.value=test.out$p.value,
			   mean._control=mean(xdata.v[g!=group,][["y"]]),
			   mean._case=mean(xdata.v[g==group,][["y"]]),
			   sd._control=sd(xdata.v[g!=group,][["y"]]),
			   sd._case=sd(xdata.v[g==group,][["y"]])
			   )
	    }
	},.parallel=T))
	all.table$adj.P.Value <- p.adjust(all.table$p.value,method="BH")
	all.table$wilcox.adj.p.value <- p.adjust(all.table$wilcox.p.value,method="BH")
	
	all.table <- all.table[order(adj.P.Value,-t,-logFC),]
	sig.table <- all.table[adj.P.Value<T.fdr & abs(logFC)>T.logFC,]

    }

    ret.dat <- list(all=all.table,sig=sig.table)
    return(ret.dat)

}


