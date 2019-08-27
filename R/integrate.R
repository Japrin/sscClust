#' integrate dataset by calculating the average expression and then clustering the pseudo bulk data
#' @param sce.list list; list of object of \code{singleCellExperiment} class
#' @param out.prefix character; output prefix.
#' @param assay.name character; which assay (default: "exprs")
#' @param is.avg logical; whether the objects in sce.list are average expression (default: FALSE)
#' @param ncores integer; number of cores to use (default: 6)
#' @param use.deg logical; whether use only the differentially expressed genes (default: FALSE)
#' @param gene.de.list list; if not NULL, each element is a data.frame with "geneID" column (default: NULL)
#' @param avg.by character; calculate the average expression of cells group by the specifid column (default: "majorCluster")
#' @param n.downsample integer; number of cells in each cluster to downsample to (default: NULL)
#' @param n.pc integer; number of pc ot use (default: 15)
#' @param method.avg character; method of calculate the average expression. Passed to `avg` of `ssc.average.cell`.(default: "zscore")
#' @param topGene.lo double; for top gene heatmap.(default: -1.5)
#' @param topGene.hi double; for top gene heatmap.(default: 1.5)
#' @param topGene.step double; for top gene heatmap.(default: 1)
#' @param ... parameters passed to ssc.clusterMarkerGene
#' @importFrom plyr llply
#' @details method to calculate the average expression can be one of "mean", "zscore"
#' @export
integrate.by.avg <- function(sce.list,
                             out.prefix,
                             assay.name="exprs",is.avg=FALSE,ncores=6,
                             use.deg=TRUE,
                             gene.de.list=NULL,
                             avg.by="majorCluster",
							 n.downsample=NULL,
							 n.pc=15,
							 ###par.clust=list(deepSplit=4, minClusterSize=2,method="dynamicTreeCut"),
							 par.clust=list(method="SNN",SNN.k=3,SNN.method="leiden",resolution_parameter=2.2),
                             topGene.lo=-1.5,topGene.hi=1.5,topGene.step=1,
                             method.avg="zscore",...)
  {
    require("plyr")
    if(use.deg && is.null(gene.de.list) && is.avg){
        cat(sprintf("sce.list contain average expression, gene.de.list must be not NULL!\n"))
    }
    if(is.null(names(sce.list))){
        names(sce.list) <- sprintf("D%02d",seq_along(sce.list))
    }

    #### get common genes
    loginfo(sprintf("get common genes ... "))
    gene.common <- c()
    for(i in seq_along(sce.list)){
      if(length(gene.common)==0)
            gene.common <- rownames(sce.list[[i]])
        else
            gene.common <- intersect(gene.common,rownames(sce.list[[i]]))
    }
    loginfo(sprintf("total %d common genes obtain ",length(gene.common)))

    ##Differential expressed genes
    gene.de.common <- c()
    if(use.deg){
		loginfo(sprintf("use deg "))
        if(is.null(gene.de.list)){
            gene.de.list <- list()
            for(i in seq_along(sce.list)){
				loginfo(sprintf(">>> cal deg for %s ",names(sce.list)[i]))
                de.out <- ssc.clusterMarkerGene(sce.list[[i]],assay.name=assay.name,
												ncell.downsample=n.downsample,
												n.cores=ncores,group.var=avg.by,...)
                gene.de.list[[i]] <- de.out$gene.table
            }
            names(gene.de.list) <- names(sce.list)
        }
        for(i in seq_along(gene.de.list)){
            if(length(gene.de.common)==0)
                gene.de.common <- gene.de.list[[i]]$geneID
            else
                gene.de.common <- c(gene.de.common,gene.de.list[[i]]$geneID)
        }
        gene.de.common <- unique(gene.de.common)
		loginfo(sprintf("total number %d deg in the union set ",length(gene.de.common)))
        #gene.common <- intersect(gene.common,gene.de.common)
    }
	if(length(gene.de.common)>0){
		gene.de.common <- intersect(gene.common,gene.de.common)
	}else{
		gene.de.common <- gene.common
	}
	loginfo(sprintf("total number %d genes will be used ",length(gene.de.common)))

    ### cal the average
    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = ncores)
	loginfo(sprintf("cal the average expression of each cluster "))
    sce.avg.list <- llply(seq_along(sce.list),function(i){
                             aid <- names(sce.list)[i]
                             obj <- sce.list[[i]][gene.common,]
                             colData(obj)[,avg.by] <- sprintf("%s.%s",aid,colData(obj)[,avg.by])
                             obj.avg <- obj
                             if(!is.avg){
                                 obj.avg <- ssc.average.cell(obj,assay.name = assay.name,
                                                             avg=method.avg,column=avg.by,ret.type="sce")
                             }else{
                                 colnames(obj.avg) <- sprintf("%s.%s",aid,colnames(obj.avg))
                             }
						     loginfo(sprintf(">>> average expression calculation done for %s (ngenes: %d, ncells: %d)",
											 aid,nrow(obj),ncol(obj)))
                             return(obj.avg)
    },.parallel = T)
    names(sce.avg.list) <- names(sce.list)

	loginfo(sprintf("get the average expression data across datasets "))
    dat.avg.mtx <- NULL
    for(aid in names(sce.avg.list)){
        if(is.null(dat.avg.mtx)){
            dat.avg.mtx <- assay(sce.avg.list[[aid]],assay.name)
        }else{
            dat.avg.mtx <- cbind(dat.avg.mtx,assay(sce.avg.list[[aid]],assay.name))
        }
    }

    ###  clustering the pseudo bulk data
    #all(rownames(dat.avg.mtx)==rownames(sce.avg.list[[1]]))
	loginfo(sprintf("cluster the pseudo-bulk samples..."))
    sce.pb <- ssc.build(dat.avg.mtx)

    rowData(sce.pb)$gene.de.common <- rownames(sce.pb) %in% gene.de.common

    sce.pb <- ssc.reduceDim(sce.pb,method="pca",method.vgene="gene.de.common",
                            pca.npc=n.pc,seed=9997)

    #ssc.plot.pca(sce.pb)
    m <- regexec("^(.+?)\\.",colnames(sce.pb),perl=T)
    mm <- regmatches(colnames(sce.pb),m)
    colData(sce.pb)[,avg.by] <- colnames(sce.pb)
    colData(sce.pb)[,"dataset.id"] <- sapply(mm,"[",2)
    colData(sce.pb)

    #### clustering the clusters
	sce.pb <- do.call(ssc.clust,c(list(obj=sce.pb,method.reduction="pca", seed=9997),
								 par.clust))

	loginfo(sprintf("make some plots ..."))
    
	p <- ssc.plot.tsne(sce.pb,columns = "dataset.id",reduced.name = "pca.tsne",size=3)
    ggsave(sprintf("%s.pca.tsne.aid.pdf",out.prefix),width=5,height=4)
	
    #p <- ssc.plot.tsne(sce.pb,columns = "pca.hclust.k0",
    #               reduced.name = "pca.tsne",
    #               size=3)
    #ggsave(sprintf("%s.pca.tsne.hclust.k0.pdf",out.prefix),width=5.5,height=4)

	if("pca.dynamicTreeCut.kauto" %in% colnames(colData(sce.pb)))
	{
		p <- ssc.plot.tsne(sce.pb,columns = "pca.dynamicTreeCut.kauto",
						   reduced.name = "pca.tsne",
						   size=3)
		ggsave(sprintf("%s.pca.tsne.dynamicTreeCut.kauto.pdf",out.prefix),
				   width=6.0,height=4)

		branch.out <- sscClust:::plot.branch(metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$hclust,
					cluster=sce.pb$pca.dynamicTreeCut.kauto,
					out.prefix=sprintf("%s.dynamicTreeCut",out.prefix))
		metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch <- branch.out$branch

		#### correlation between clusters
		dat.cor.mtx <- cor(assay(sce.pb),method="pearson")
		for(z.lo in c(0,0.5,-0.5,-1)){
			sscClust:::plot.matrix.simple(dat.cor.mtx,
							   out.prefix=sprintf("%s.avg.%s.cor.zlo.%s",out.prefix,"pearson",z.lo),
							   do.clust=metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch,
							   show.number = F,
							   z.lo = z.lo, z.hi = 1,exp.name="Cor")
		}

		sscClust:::plot.matrix.simple(as.matrix(metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$dist),
						   out.prefix=sprintf("%s.avg.pca.dist",out.prefix),
						   do.clust=metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch,show.number = F,
						   palatte=(brewer.pal(n = 7,name = "RdYlBu")),
						   z.lo = 0, z.hi = 15,exp.name="dist")
	}

	if(!is.null(gene.de.list) && "pca.SNN.kauto" %in% colnames(colData(sce.pb))){
		##### top de genes
		gene.desc.top <- as.data.table(ldply(names(gene.de.list),function(aid){
								   gene.de.list[[aid]]$Group <- sprintf("%s.%s",aid,gene.de.list[[aid]][["cluster"]])
								   return(gene.de.list[[aid]][,c("geneID","geneSymbol","t","cluster","Group"),with=F])
								   }))
		gene.desc.top$meta.cluster <- sce.pb$pca.SNN.kauto[match(gene.desc.top$Group,colnames(sce.pb))]
		gene.desc.top <- gene.desc.top[order(meta.cluster,-t,Group),]
		saveRDS(gene.desc.top,sprintf("%s.gene.desc.top.rds",out.prefix))
		##gene.desc.top <- readRDS(sprintf("%s.gene.desc.top.rds",out.prefix))

		RhpcBLASctl::omp_set_num_threads(1)
		doParallel::registerDoParallel(cores = ncores)

		g.desc <- ldply(sort(unique(gene.desc.top$meta.cluster)),function(mcls){
					dat.plot <- gene.desc.top[meta.cluster==mcls,]
					ncluster <- length(unique(sort(dat.plot$Group)))
					gene.core.tb <- dat.plot[,.(N=.N),by=c("meta.cluster","geneID")][N>0.3*ncluster,]
					gene.core.tb$meta.cluster.size <- ncluster
					gene.core.tb$Group <- gene.core.tb$meta.cluster
					g.desc <- head(gene.core.tb[geneID %in% rownames(sce.pb),],n=30)
					ssc.plot.heatmap(sce.pb,out.prefix=sprintf("%s.gene.top.meta.cluster.%s",out.prefix,mcls),
							   columns="pca.SNN.kauto",columns.order="pca.SNN.kauto",
							   gene.desc=g.desc,
							   pdf.width=20,pdf.height=10,do.scale=F,
							   z.lo=topGene.lo,z.hi=topGene.hi,z.step=topGene.step,
							   do.clustering.row=F,
							   do.clustering.col=T
							   )
					return(g.desc)
				   },.parallel=T)

		ssc.plot.heatmap(sce.pb,out.prefix=sprintf("%s.gene.top.meta.cluster.all",out.prefix),
							   columns="pca.SNN.kauto",columns.order="pca.SNN.kauto",
							   gene.desc=as.data.table(g.desc)[!duplicated(geneID),],
							   pdf.width=20,pdf.height=15,do.scale=F,
							   z.lo=topGene.lo,z.hi=topGene.hi,z.step=topGene.step,
							   do.clustering.row=F,
							   do.clustering.col=T
							   )
	}
	
	#### heatmap show average expression of specified genes
    ###return(list("sce.pb"=sce.pb,"branch.out"=branch.out))
    metadata(sce.pb)$ssc$gene.de.list <- gene.de.list
    metadata(sce.pb)$ssc$gene.desc.top <- gene.desc.top
	loginfo(sprintf("integrate.by.avg run successfully"))
    return(sce.pb)
}

#' plot genes expression in pairs of clusters to examine the correlation
#' @param obj.list object; named list of object of \code{singleCellExperiment} class
#' @param gene.desc.top data.frame; signature genes 
#' @param assay.name character; which assay (default: "exprs")
#' @param out.prefix character; output prefix (default: NULL).
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param sig.prevelance double; (default: 0.5)
#' @import data.table
#' @import ggplot2 
#' @import ggridges
#' @importFrom  ggpubr ggscatter
#' @details classify cells using the signature genes
#' @export
classifyCell.by.sigGene <- function(obj.list,gene.desc.top,assay.name="exprs",out.prefix=NULL,
                                    adjB=NULL,
                                    sig.prevelance=0.5){
    mcls <- unique(gene.desc.top$meta.cluster)

    dat.list <- (llply(seq_along(mcls),function(j){
        gene.desc.top.j <- gene.desc.top[meta.cluster==mcls[j],]
		ncluster <- length(unique(sort(gene.desc.top.j$Group)))
		gene.core.tb <- gene.desc.top.j[,.(N=.N),by=c("meta.cluster","geneID")][N>sig.prevelance*ncluster,]
		gene.core.tb$meta.cluster.size <- ncluster
		gene.core.tb$Group <- gene.core.tb$meta.cluster
		dat.plot.j <- ldply(seq_along(obj.list),function(i){
							  obj <- obj.list[[i]]
							  gene.used <- rownames(obj)[match(gene.core.tb$geneID,rowData(obj)$display.name)]
							  gene.used <- gene.used[!is.na(gene.used)]
							  dat.block <- assay(obj,assay.name)[gene.used,,drop=F]
							  if(!is.null(adjB)){
								dat.block <- simple.removeBatchEffect(dat.block,batch=obj[[adjB]])
							  }
							  sig.score <- colMeans(dat.block)
							  data.table(cellID=colnames(obj),
										 dataset.id=names(obj.list)[i],
										 meta.cluster=mcls[j],
										 sig.score=sig.score)
									})
		return(list(dat.plot=dat.plot.j,gene.core.tb=gene.core.tb))
	}))
	dat.plot <- as.data.table(ldply(dat.list,function(x){ x$dat.plot }))
	gene.core.tb <- as.data.table(ldply(dat.list,function(x){ x$gene.core.tb }))

    if(!is.null(out.prefix)){
		p <- ggplot(dat.plot, aes(x=sig.score,y=dataset.id,color=dataset.id),fill="none") +
			##geom_histogram(aes(y=..density..), colour="black", fill="white")+
			#geom_density() +
			geom_density_ridges() +
			facet_wrap(~meta.cluster,ncol=4)
		ggsave(file=sprintf("%s.sigScore.density.pdf",out.prefix),width=15,height=12)
    }

    
}

####
#####' outlier detection using extremevalues
#####' @param x object; vector
#####' @import data.table
#####' @import ggplot2 
#####' @import ggridges
#####' @import extremevalues
#####' @importFrom  ggpubr ggscatter
#####' @details outlier detection using extremevalues
#####' @export
####classify.outlier <- function(x){
####	##### outlier detection (use method I at last)
####	K <- getOutliers(x,method="I",distribution="normal")
####	L <- getOutliers(x,method="II",distribution="normal")
####	pdf(sprintf("%s.outlier.pdf",out.prefix),width = 10,height = 6)
####	opar <- par(mfrow=c(1,2))
####	outlierPlot(myDesign$prol.score,K,mode="qq")
####	outlierPlot(myDesign$prol.score,L,mode="residual")
####	dev.off()
####	par(opar)
####	print("K$mu,K$sigma,K$limit")
####	print(c(K$mu,K$sigma,K$limit))
####	myDesign$prol.cls <- "lo"
####	myDesign$prol.cls[K$iRight] <- "hi"
####	######
####
####	##### plot #####
####	#prol.score.x <- seq(-2,6,0.001)
####	##prol.score.x <- seq(K$yMin,K$yMax,0.001)
####	prol.score.pretty <- pretty(myDesign$prol.score)
####	prol.score.x <- seq(prol.score.pretty[1], prol.score.pretty[length(prol.score.pretty)],0.001)
####	prol.score.forPlot <- data.frame(x=prol.score.x,y=dnorm(prol.score.x,mean = K$mu,sd = K$sigma))
####	p <- ggplot(myDesign, aes(prol.score)) + geom_density(colour="black") + theme_bw() +
####		geom_line(data=prol.score.forPlot,aes(x=x,y=y),colour="red") +
####		geom_vline(xintercept = K$limit,linetype=2)
####	ggsave(filename = sprintf("%s.density.pdf",out.prefix),width = 4,height = 3)
####
####
####}
####
####### zero.as.low: if True, zero counts as low expression (0); else as "-1"
####binarizedExp <- function(x,ofile=NULL,G=NULL,e.TH=NULL,e.name="Exp",verbose=F, draw.CI=T, zero.as.low=T,...)
####{
####  require(mclust)
####  ### bin.Exp == -1: drop out
####  o.df <- data.frame(sample=names(x),bin.Exp=-1,stringsAsFactors = F)
####  rownames(o.df) <- o.df$sample
####  #f<-is.finite(x) & x>0
####  f<-is.finite(x)
####  if(sum(f)<3){
####      colnames(o.df) <- c("sample",e.name)
####      if(verbose){
####        return(list(x_mix=NULL,o.df=o.df))
####      }else{
####        return(structure(o.df[,e.name],names=rownames(o.df)))
####      }
####  }
####
####  x<-x[f]
####  x_mix<-densityMclust(x,G=G,modelNames=c("E","V"))
####  x_mix_summary<-summary(x_mix)
####  
####  if(verbose && !is.null(ofile)){
####	  print(x_mix_summary)
####      if(grepl(".png$",ofile,perl = T)){
####          png(ofile,width=800,height=600)
####	      old_par<-par(no.readonly=T)
####	      layout(matrix(c(1,1,1,1,1,1,2,3,4), 3, 3, byrow = TRUE))
####	      a_par<-par(cex.axis=2,cex.lab=2,cex.main=1.8,mar=c(5,6,4,2)+0.1)
####      }else{
####          pdf(ofile,width=8,height=6)
####      }
####	  plot(x_mix,what="density",data=x,breaks=50,col="darkgreen",lwd=2,main="",...)
####	  abline(v=x_mix_summary$mean,lty=2)
####	  if(!is.null(e.TH)){
####		abline(v=e.TH,lty=2,col="red")
####	  }
####	  
####	  for(i in 1:x_mix_summary$G)
####	  {
####		i_mean<-x_mix_summary$mean[i]
####		i_sd<-sqrt(x_mix_summary$variance[i])
####		i_pro<-x_mix_summary$pro[i]
####		#i_sd<-RC_mix_summary$variance[i]
####		d<-qnorm(c(0.0013,0.9987),i_mean,i_sd)
####		e<-i_pro*dnorm(i_mean,i_mean,i_sd)
####		lines(seq(d[1],d[2],by=0.01),i_pro*dnorm(seq(d[1],d[2],by=0.01),i_mean,i_sd),col="orange",lwd=2)
####        if(draw.CI){ rect(d[1],0,d[2],e+0.02,col=NA,border="blue") }
####		#rect(d[1],0,d[2],e+0.02,col=rgb(0,0,0.8,0.2),border=NA)
####	  }
####	  plot(x_mix,data=x,breaks=20,col="darkgreen",lwd=2,what="BIC")
####	  densityMclust.diagnostic(x_mix,type = "cdf",cex.lab=1.5)
####	  densityMclust.diagnostic(x_mix,type = "qq")
####	  dev.off()
####  }
####  
####  o.df[names(x_mix$classification),"bin.Exp"] <- x_mix$classification 
####  colnames(o.df) <- c("sample",e.name)
####  if(!is.null(G) && x_mix_summary$G==3){
####	  ### determin which classes 'not-expressed' and which classes 'expressed'	
####	  i_mean<-x_mix_summary$mean
####	  i_sd<-sqrt(x_mix_summary$variance)
####	  ci95.1 <- qnorm(c(0.0013,0.9987),i_mean[1],i_sd[1])
####	  ci95.1.overlap <- pnorm(ci95.1[2],i_mean[2],i_sd[2])
####	  if(ci95.1.overlap>0.3333){
####		  C.min <-  2
####	  }else{
####		  C.min <- 1
####	  }
####	  ci95.2 <- qnorm(c(0.0013,0.9987),i_mean[2],i_sd[2])
####	  ci95.3 <- qnorm(c(0.0013,0.9987),i_mean[3],i_sd[3])
####	  ci95.3.overlap <- pnorm(ci95.3[1],i_mean[2],i_sd[2],lower.tail = F)
####	  if(ci95.3.overlap>0.3333){
####		  C.max <- 2
####	  }else{
####		  C.max <- 3
####	  }
#####	  ci95.2.overlap <- pnorm(ci95.2[2],i_mean[3],i_sd[3])
#####	  if(ci95.2.overlap>0.3333){
#####	  	  C.max <- 2
#####	  }else{
#####	  	  C.max <- 3
#####	  }
####      if(zero.as.low){
####	    f.low <- o.df[,e.name] <= C.min
####      }else{
####	    f.low <- o.df[,e.name] <= C.min & o.df[,e.name] != -1
####      }
####	  f.hi <-  o.df[,e.name] >= C.max
####	  f.mid <- (o.df[,e.name] > C.min) & (o.df[,e.name] < C.max)
####	  o.df[f.low,e.name] <- 0
####	  o.df[f.mid,e.name] <- 0.5
####	  o.df[f.hi,e.name] <- 1
####	  #f.G <- o.df[,e.name] < C.max
####	  #o.df[f.G,e.name] <- 0
####	  #o.df[!f.G,e.name] <- 1
####  }else if(x_mix_summary$G==2){
####      if(zero.as.low){
####	    f.low <- o.df[,e.name] <= 1
####      }else{
####	    f.low <- o.df[,e.name] <= 1 & o.df[,e.name] != -1
####      }
####	  f.hi <-  o.df[,e.name] == 2
####	  o.df[f.low,e.name] <- 0
####	  o.df[f.hi,e.name] <- 1
####  }
####  if(verbose){
####	return(list(x_mix=x_mix,o.df=o.df))
####  }else{
####	return(structure(o.df[,e.name],names=rownames(o.df)))
####  }
####}
####

#' plot genes expression in pairs of clusters to examine the correlation
#' @param sce.pb object; object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param x.cluster character; cluster to be showed in x axis
#' @param out.prefix character; output prefix.
#' @param genes vector; genes to used to calculate the correlation and make the scatter plot (default: NULL)
#' @param gene.highlight vector; genes to highlight (default: NULL)
#' @param plot.ncol integer; number of columns in facet_wrap() (default: 7)
#' @param plot.width double; width of the output plot (default: 22)
#' @param plot.height double; height of the output plot (default: 22)
#' @import data.table
#' @import ggplot2
#' @importFrom  ggpubr ggscatter
#' @details used to examine the correlation between pairs of clusters
#' @export
plot.pairs.cor <- function(sce.pb,assay.name,x.cluster,out.prefix,
                           genes=NULL,gene.highlight=NULL,
                           plot.ncol=7,plot.width=22,plot.height=22)
{
    require("data.table")
    require("ggplot2")
    require("ggpubr")
    dat.assay <- assay(sce.pb,assay.name)
    colnames(dat.assay) <- make.names(colnames(dat.assay))
    if(!is.null(genes)){
        f.gene <- intersect(rownames(dat.assay),genes)
        dat.assay <- dat.assay[f.gene,]
    }
    dat.plot <- setDT(as.data.frame(dat.assay),keep.rownames=T)
    dat.plot <- melt(dat.plot,id=c("rn",x.cluster))
    dat.plot$isMarker <- FALSE
    if(!is.null(gene.highlight)){
        dat.plot$isMarker <- dat.plot$rn %in% gene.highlight
    }
    dat.plot[1:4,]

    dat.cor <- cor(dat.assay)
    set.cls <- sort(dat.cor[,x.cluster],decreasing = T)
    dat.plot$variable <- factor(as.character(dat.plot$variable), levels=names(set.cls))

    p <- ggscatter(dat.plot,x=x.cluster,y="value",
                   fill="isMarker",color="isMarker",palette=c("FALSE"="lightblue","TRUE"="red"),
                   add = "reg.line",
                   add.params = list(color = "blue", fill = "lightgray"),
                   conf.int = TRUE,
                   cor.coef = TRUE,
                   cor.coeff.args = list(method = "pearson", label.x = -0.5, label.sep = "\n"),
                   size=c(0.5,2)[as.integer(dat.plot$isMarker)+1]) +
            geom_hline(yintercept=c(-1,-0.5,0.5,1),linetype="dashed",color="gray") +
            geom_vline(xintercept=c(-1,-0.5,0.5,1),linetype="dashed",color="gray") +
            facet_wrap(~variable,ncol=plot.ncol,scales="free")
    ggsave(sprintf("%s.%s.png",out.prefix,x.cluster),width=plot.width,height=plot.height)
}


