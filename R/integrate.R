#' integrate dataset by calculating the average expression and then clustering the pseudo bulk data
#' @param sce.list list; list of object of \code{singleCellExperiment} class
#' @param out.prefix character; output prefix.
#' @param assay.name character; which assay (default: "exprs")
#' @param is.avg logical; whether the objects in sce.list are average expression (default: FALSE)
#' @param ncores integer; number of cores to use (default: 6)
#' @param use.deg logical; whether use only the differentially expressed genes (default: FALSE)
#' @param gene.de.list list; if not NULL, each element is a data.frame with "geneID" column (default: NULL)
#' @param sort.by character; by which to sort genes and the top ones will be used for clustering. One of "F.rank.median", and "occurence" (default: "F.rank.median")
#' @param avg.by character; calculate the average expression of cells group by the specifid column (default: "majorCluster")
#' @param n.downsample integer; number of cells in each cluster to downsample to (default: NULL)
#' @param n.pc integer; number of pc ot use (default: 15)
#' @param do.clustering logical; wheter perform PCA/clustering (default: TRUE)
#' @param de.stat character; column in gene.de.file (default: "t")
#' @param de.thres double; used for selecting genes. can be larger than 1(top number of genes) or less than 1 (F.rank threshold) (default: 2000)
#' @param do.scale logical; scale the summarized vectors (default: false)
#' @param par.clust list; parameters for clustering method (default: list(method="SNN",SNN.k=3,SNN.method="leiden",resolution_parameter=2.2))
#' @param method.avg character; method of calculate the average expression. Passed to `avg` of `ssc.average.cell`.(default: "zscore")
#' @param topGene.lo double; for top gene heatmap.(default: -1.5)
#' @param topGene.hi double; for top gene heatmap.(default: 1.5)
#' @param topGene.step double; for top gene heatmap.(default: 1)
#' @param myseed integer; seed for random number generation.(default: 9997)
#' @param verbose logical; verbose level.(default: FALSE)
#' @param ... parameters passed to ssc.clusterMarkerGene
#' @importFrom plyr llply
#' @importFrom SummarizedExperiment colData rowData `colData<-` `rowData<-` assay
#' @importFrom S4Vectors `metadata<-`
#' @importFrom matrixStats rowMedians
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils head
#' @details method to calculate the average expression can be one of "mean", "zscore"
#' @export
integrate.by.avg <- function(sce.list,
                             out.prefix,
                             assay.name="exprs",is.avg=FALSE,ncores=6,
                             use.deg=TRUE,
                             gene.de.list=NULL,sort.by="F.rank.median",
                             avg.by="majorCluster",
							 n.downsample=NULL,
							 n.pc=15,do.clustering=T,
							 de.stat="t",de.thres=1,do.scale=F,
							 ###par.clust=list(deepSplit=4, minClusterSize=2,method="dynamicTreeCut"),
							 par.clust=list(method="SNN",SNN.k=3,SNN.method="leiden",resolution_parameter=2.2),
                             topGene.lo=-1.5,topGene.hi=1.5,topGene.step=1,myseed=9997,verbose=F,
                             method.avg="zscore",...)
  {
    #require("plyr")
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
	sce.pb <- NULL
	for(xx in assayNames(sce.avg.list[[1]])){
		dat.avg.mtx <- NULL
		for(aid in names(sce.avg.list)){
			if(is.null(dat.avg.mtx)){
				dat.avg.mtx <- assay(sce.avg.list[[aid]],xx)
			}else{
				dat.avg.mtx <- cbind(dat.avg.mtx,assay(sce.avg.list[[aid]],xx))
			}
		}
		if(is.null(sce.pb)){
			sce.pb <- ssc.build(dat.avg.mtx,assay.name=xx)
		}else{
			assay(sce.pb,xx) <- dat.avg.mtx
		}
	}
	assay(sce.pb,"exprs") <- assay(sce.pb,assay.name)
	##if(is.avg)
	{
		cls.info.df <- ldply(names(sce.avg.list),function(x){
								 o.df <- cbind(data.frame(rid=colnames(sce.avg.list[[x]])),
											   as.data.frame(colData(sce.avg.list[[x]])))
								 o.df
								})
		rownames(cls.info.df) <- cls.info.df$rid
		colData(sce.pb) <- DataFrame(cls.info.df[colnames(sce.pb),])
	}

    ##Differential expressed genes
    gene.de.common <- c()
	#use.F.avg.rank <- TRUE
    if(use.deg && !is.null(gene.de.list)){
		if(sort.by=="F.rank.median"){
			gene.rank.tb <- as.data.table(ldply(names(gene.de.list),function(x){
													ret.tb <- unique(gene.de.list[[x]][,c("geneID","F.rank")])
													ret.tb$dataset.id <- x
													return(ret.tb)
															 }))
			gene.rank.tb <- dcast(gene.rank.tb,geneID~dataset.id,value.var="F.rank",fill=1)
			gene.rank.tb$median.F.rank <- rowMedians(as.matrix(gene.rank.tb[,-c("geneID"),with=F]),na.rm=T)
			gene.rank.tb <- gene.rank.tb[order(median.F.rank),]
			gene.rank.tb <- gene.rank.tb[geneID %in% gene.common,]
			rowData(sce.pb)$median.F.rank <- gene.rank.tb[["median.F.rank"]][match(rownames(sce.pb),gene.rank.tb$geneID)]
			#gene.rank.tb.debug <<- gene.rank.tb
			#sce.pb.debug <<- sce.pb
			if(de.thres>=1){
				gene.de.common <- head(gene.rank.tb$geneID,n=de.thres)
			}else{
				gene.de.common <- gene.rank.tb[median.F.rank < de.thres,][["geneID"]]
			}
			if(verbose){

#				.dat.plot <- gene.rank.tb
#				.dat.plot$rank <- order(gene.rank.tb$median.F.rank)
#				.dat.plot$used <- FALSE
#				.dat.plot[rank<=de.thres,used:=TRUE]
#				p <- ggplot(.dat.plot,aes(x=rank,y=median.F.rank)) +
#						geom_point(aes(color=used),shape=20) +
#						xlab("Genes") + ylab("Median of\nPercentage Rank") +
#						theme_bw() + 
#						coord_flip() + scale_x_reverse() +
#						theme(axis.title=element_text(size=18)) + 
#						theme(axis.text=element_text(size=15))
#				if(de.thres>=1){ p <- p + geom_vline(xintercept=de.thres,linetype=2,color="black") }
#				#ggsave(sprintf("%s.Frank.rank.png",out.prefix),width=6,height=2.5)
#				ggsave(sprintf("%s.Frank.rank.png",out.prefix),width=4,height=6)

			}
		}else if(sort.by=="occurence")
		{
			loginfo(sprintf("use deg present multiple times (de.thres: %s) ",de.thres))
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
					gene.de.common <- unique(gene.de.list[[i]]$geneID)
				else
					gene.de.common <- c(gene.de.common,unique(gene.de.list[[i]]$geneID))
			}
			#gene.de.common <- unique(gene.de.common)
			de.gene.ntimes <- table(gene.de.common)
			if(de.thres>=1){
				gene.de.common <- names(de.gene.ntimes[de.gene.ntimes >= de.thres])
			}else{
				gene.de.common <- names(de.gene.ntimes[de.gene.ntimes >= de.thres*length(gene.de.list)])
			}
			loginfo(sprintf("total number %d deg in the union set ",length(gene.de.common)))
		}
    }
	if(length(gene.de.common)>0){
		gene.de.common <- intersect(gene.common,gene.de.common)
	}else{
		gene.de.common <- gene.common
	}
	loginfo(sprintf("total number %d genes will be used ",length(gene.de.common)))

	rowData(sce.pb)$gene.de.common <- rownames(sce.pb) %in% gene.de.common

    m <- regexec("^(.+?)\\.",colnames(sce.pb),perl=T)
    mm <- regmatches(colnames(sce.pb),m)
    colData(sce.pb)[,avg.by] <- colnames(sce.pb)
    colData(sce.pb)[,"dataset.id"] <- sapply(mm,"[",2)
    colData(sce.pb)
    metadata(sce.pb)$ssc$gene.de.list <- gene.de.list

	if(!do.clustering){
		loginfo(sprintf("integrate.by.avg run successfully"))
		return(sce.pb)
	}

    ###  clustering the pseudo bulk data
    #all(rownames(dat.avg.mtx)==rownames(sce.avg.list[[1]]))
	loginfo(sprintf("cluster the pseudo-bulk samples..."))

	if(do.scale){
		assay(sce.pb,"exprs.scale") <- t(scale(t(assay(sce.pb,"exprs"))))
		sce.pb <- ssc.reduceDim(sce.pb,assay.name="exprs.scale",
								method="pca",method.vgene="gene.de.common", pca.npc=n.pc,seed=myseed)
	}else{
		sce.pb <- ssc.reduceDim(sce.pb,assay.name="exprs",
								method="pca",method.vgene="gene.de.common", pca.npc=n.pc,seed=myseed)
	}

    #ssc.plot.pca(sce.pb)
    #### clustering the clusters
	sce.pb <- do.call(ssc.clust,c(list(obj=sce.pb,method.reduction="pca", seed=myseed),
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

		branch.out <- sscVis:::plotBranch(metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$hclust,
					cluster=sce.pb$pca.dynamicTreeCut.kauto,
					out.prefix=sprintf("%s.dynamicTreeCut",out.prefix))
		metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch <- branch.out$branch

		#### correlation between clusters
		dat.cor.mtx <- cor(assay(sce.pb),method="pearson")
		for(z.lo in c(0,0.5,-0.5,-1)){
			sscVis::plotMatrix.simple(dat.cor.mtx,
							   out.prefix=sprintf("%s.avg.%s.cor.zlo.%s",out.prefix,"pearson",z.lo),
							   do.clust=metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch,
							   show.number = F,
							   z.lo = z.lo, z.hi = 1,exp.name="Cor")
		}

		sscVis::plotMatrix.simple(as.matrix(metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$dist),
						   out.prefix=sprintf("%s.avg.pca.dist",out.prefix),
						   do.clust=metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch,show.number = F,
						   palatte=(brewer.pal(n = 7,name = "RdYlBu")),
						   z.lo = 0, z.hi = 15,exp.name="dist")
	}

	if(!is.null(gene.de.list) && "pca.SNN.kauto" %in% colnames(colData(sce.pb))){
		##### top de genes
		gene.desc.top <- rank.de.gene(sce.pb)
		metadata(sce.pb)$ssc$gene.desc.top <- gene.desc.top
		#saveRDS(gene.desc.top,sprintf("%s.gene.desc.top.rds",out.prefix))
		##gene.desc.top <- readRDS(sprintf("%s.gene.desc.top.rds",out.prefix))

		RhpcBLASctl::omp_set_num_threads(1)
		doParallel::registerDoParallel(cores = ncores)

		g.desc <- ldply(sort(unique(gene.desc.top$meta.cluster)),function(mcls){
					dat.plot <- gene.desc.top[meta.cluster==mcls,]
					gene.core.tb <- dat.plot[freq.sig>0.3 & median.rank < 0.1,]
#					ncluster <- length(unique(sort(dat.plot$Group)))
#					gene.core.tb <- dat.plot[,.(N=.N),by=c("meta.cluster","geneID")][N>0.3*ncluster,]
#					gene.core.tb$meta.cluster.size <- ncluster
#					gene.core.tb$Group <- gene.core.tb$meta.cluster
					g.desc <- head(gene.core.tb[geneID %in% rownames(sce.pb),],n=30)
					ssc.plot.heatmap(sce.pb,out.prefix=sprintf("%s.gene.top.meta.cluster.%s",out.prefix,mcls),
							   columns="pca.SNN.kauto",columns.order="pca.SNN.kauto",
							   gene.desc=g.desc,mytitle=sprintf("%s",mcls),
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
							   mytitle=sprintf("all meta-clusters"),
							   z.lo=topGene.lo,z.hi=topGene.hi,z.step=topGene.step,
							   do.clustering.row=F,
							   do.clustering.col=T
							   )
	}

	loginfo(sprintf("integrate.by.avg run successfully"))
    return(sce.pb)
}


#' collapse effect size, give input table and group variate
#' @param dat.long data.table; these columns are required: "geneID","aid","P.Value","adj.P.Val","sig","dprime","vardprime","group.2nd"
#' @param group.2nd character; column of obj, used for collapse [default: NULL]
#' @param mode.collapse character; collapse method. one of "min", "comb" [default: "comb"]
#' @param th.adj.P double; threshold [default: 0.01]
#' @param th.dprime double; threshold [default: 0.15]
#' @importFrom data.table as.data.table data.table setorderv
#' @importFrom plyr ldply
#' @details collapse effect size
#' @return a gene table
collapseEffectSizeLong <- function(dat.long,mode.collapse="comb",group.2nd="group.2nd",
								   th.adj.P=0.01,th.dprime=0.15)
{
	#dat.long[["group.2nd"]] <- dat.long[[group.2nd]]
	dat.collapse <- NULL
	if(mode.collapse=="min"){
		dat.collapse <- dat.long[,.(dprime=dprime[which.min(abs(dprime))],
									vardprime=vardprime[which.min(abs(dprime))],
									P.Value=P.Value[which.min(abs(dprime))],
									adj.P.Val=adj.P.Val[which.min(abs(dprime))],
									sig=all(adj.P.Val<th.adj.P) & all(dprime>th.dprime) ),
								by=c("geneID",group.2nd)]
	}else if(mode.collapse=="comb"){
		aid.mapping.tb <- unique(dat.long[,c("aid",group.2nd),with=F])
		n.group.2nd <- aid.mapping.tb[,.N,by=group.2nd]
		dat.collapse.nOne <- dat.long[dat.long[[group.2nd]] %in% n.group.2nd[N==1,][[group.2nd]],
									  c("geneID",group.2nd,"dprime","vardprime","P.Value","adj.P.Val","sig"),with=F]
		dat.collapse.nMul <- as.data.table(ldply(n.group.2nd[N>1,][[group.2nd]],
												 function(x){
													 dat.x.es.combi.tb <- directEScombiFromLongTable(dat.long[dat.long[[group.2nd]]==x,])
													 ret.tb <-data.table(geneID=dat.x.es.combi.tb$geneID,
															   group.2nd=x,
															   dprime=dat.x.es.combi.tb[["comb.ES"]],
															   vardprime=dat.x.es.combi.tb[["comb.ES.sd"]]^2,
															   P.Value=dat.x.es.combi.tb[["comb.p"]],
															   adj.P.Val=dat.x.es.combi.tb[["comb.padj"]],
															   sig=(dat.x.es.combi.tb[["comb.ES"]]>th.dprime &
																	dat.x.es.combi.tb[["comb.padj"]]<th.adj.P))
													 colnames(ret.tb)[2] <- group.2nd
													 return(ret.tb)
												 }))
		dat.collapse <- rbind(dat.collapse.nOne,dat.collapse.nMul)
		setorderv(dat.collapse,c("geneID",group.2nd))
	}
	return(dat.collapse)
}


#' get gene ranking table from obj
#' @param obj object; object of \code{singleCellExperiment} class
#' @param group character; columan name in colData(obj)
#' @param sort.by character; by which to sort the genes [default: "median.rank"]
#' @param weight.adj character; column of obj, used to adjust weight [default: NULL]
#' @param group.2nd character; column of obj, used for collapse [default: NULL]
#' @param mode.collapse character; collapse method. one of "min", "comb" [default: "comb"]
#' @param th.adj.P double; threshold [default: 0.01]
#' @param th.dprime double; threshold [default: 0.15]
#' @importFrom data.table as.data.table data.table setkey
#' @importFrom plyr ldply
#' @importFrom SummarizedExperiment assayNames
#' @importFrom stats pnorm p.adjust
#' @importFrom matrixStats rowMedians
#' @details rank genes
#' @return a gene table
rank.de.gene <- function(obj,group="pca.SNN.kauto",sort.by="median.rank",weight.adj=NULL,
						 group.2nd=NULL,mode.collapse="comb",th.adj.P=0.01,th.dprime=0.15)
{
	#### To do: add diff between this - max_other(not this)
	gene.desc.top <- as.data.table(ldply(unique(sort(obj[[group]])),function(x){
											 obj.x <- obj[,obj[[group]]==x]

											 dat.long <- ssc.toLongTable(obj.x,gene.id=NULL,
																		 assay.name=c("dprime","vardprime",
																					  "P.Value","adj.P.Val","sig"),
																		 col.idx=group.2nd)

#											 if(!is.null(group.2nd)){
#												 dat.long$group.2nd <- dat.long[[group.2nd]]
#											 }else{
#												 dat.long[,group.2nd:="GRP"]
#											 }

											 dat.collapse <- NULL

											 ret.tb <- data.table(geneID=rownames(obj.x),
														geneSymbol=rowData(obj.x)$display.name,
														meta.cluster=x,
														Group=x)

											 ret.tb.ext <- dat.long[,.(N.sig=sum(.SD$sig),
																	   freq.sig=sum(.SD$sig)/nrow(.SD),
																	   N.total=.N),
																	by="geneID"]

											 if("dprime" %in% assayNames(obj.x))
											 {
												 ret.tb.ext.a <- dat.long[,.(dprime.max=max(dprime),
																	  dprime.max.P.Value=P.Value[which.max(dprime)],
																	  dprime.max.adj.P.Val=adj.P.Val[which.max(dprime)]),
																	by="geneID"]
												
												 if(!is.null(group.2nd)){
													 dat.collapse <- collapseEffectSizeLong(dat.long,mode.collapse=mode.collapse,
																			group.2nd=group.2nd,
																			th.adj.P=th.adj.P,th.dprime=th.dprime)
												 }

												 if(!is.null(dat.collapse)){
													ret.tb.ext.b <- dat.collapse[,.(N.sig.group.2nd=sum(.SD$sig),
																			   freq.sig.group.2nd=sum(.SD$sig)/nrow(.SD),
																			   N.total.group.2nd=.N,
																			   dprime.max.group.2nd=max(dprime),
																			   dprime.max.P.Value.group.2nd=P.Value[which.max(dprime)],
																			   dprime.max.adj.P.group.2nd=adj.P.Val[which.max(dprime)]),
																			by="geneID"]
													ret.tb.ext.a <- merge(ret.tb.ext.a,ret.tb.ext.b,by="geneID")
												 }
												 ret.tb.ext <- merge(ret.tb.ext,ret.tb.ext.a,by="geneID")
											 }

											 ret.tb <- merge(ret.tb,ret.tb.ext,by="geneID")

											 ###anames <- intersect(assayNames(obj.x),c("logFC","t","SNR","meanExp","meanScale"))
											 anames <- intersect(assayNames(obj.x),c("logFC","meanScale"))
											 ret.m <- cbind(data.table(geneID=rownames(obj.x)),
															as.data.table(llply(anames,
																				function(aa){ rowMedians(assay(obj.x,aa)) }))
															)
											 colnames(ret.m)[-1] <- sprintf("median.%s",anames)
											 ret.m[["median.rank"]] <- rowMedians(apply(assay(obj.x,"t"),2,
																					 function(x){ rank(-x)/length(x)}))
											 ret.tb <- merge(ret.tb,ret.m,by="geneID")
											 setkey(ret.tb,"geneID")
											 ret.tb <- ret.tb[rownames(obj.x),]
											 all(ret.tb$geneID==rownames(obj.x))

											 if("zp" %in% assayNames(obj.x)){
												 w <- sqrt(obj.x$nCellsStudy/sum(obj.x$nCellsStudy))
												 if(!is.null(weight.adj)){
													w <- w * obj.x[[weight.adj]]
												 }
												 zsum <- (assay(obj.x,"zp") %*% w)[,1]
												 ### two-sided p value
												 #ret.tb[["p.comb.Stouffer.z"]] <- zsum
												 ret.tb[["p.comb.Stouffer"]] <- 2*(pnorm(-abs(zsum)))
												 ret.tb[["p.comb.Stouffer.adj"]] <- p.adjust(ret.tb[["p.comb.Stouffer"]],"BH")
											 }

											 if("dprime" %in% assayNames(obj.x)){
												es.combi <- directEScombi(assay(obj.x,"dprime"), assay(obj.x,"vardprime"))
												es.combi <- es.combi[ret.tb$geneID, ]
												if(!all(ret.tb$geneID==rownames(es.combi))){
													stop(sprintf("strange thing: !all(ret.tb$geneID==rownames(es.combi)) (%s)\n",
																 x))
												}
												ret.tb <- cbind(ret.tb,es.combi)
											 }

											 return(ret.tb)
							   }))
	gene.desc.top <- gene.desc.top[ order(gene.desc.top$meta.cluster,gene.desc.top[[sort.by]]), ]
	###gene.desc.top <- gene.desc.top[ order(gene.desc.top$meta.cluster,median.rank), ]
	return(gene.desc.top)
}


#' classify cells using signature genes
#' @param obj.list object; named list of object of \code{singleCellExperiment} class
#' @param gene.core.tb data.frame; signature genes 
#' @param pca.rotation matrix; rotation matrix (default: NULL)
#' @param gene.pos.tb data.frame; signature genes (positive)
#' @param gene.neg.tb data.frame; signature genes (negative)
#' @param meta.info.tb data.frame; cell information table (default: NULL)
#' @param meta.train.tb data.frame; cell information table (default: NULL)
#' @param ndownsample integer; number of cells (default: 1500)
#' @param myseed integer; seed for random number generation.(default: 123456)
#' @param assay.name character; which assay (default: "exprs")
#' @param out.prefix character; output prefix (default: NULL).
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param meta.cluster character; (default: NULL)
#' @param method character; (default: "posFreq")
#' @param verbose logical; (default: FALSE)
#' @param sig.prevelance double; (default: 0.5)
#' @param ntop integer; only use the ntop top genes (default: 10)
#' @param ncores integer; number of CPU to used (default: 16)
#' @param bin.z.pos.th numeric; z score of 'epressed' genes must be larger than this value  (default: 0.3)
#' @param bin.z.neg.th numeric; z score of 'not epressed' genes must be less than this value  (default: 0.3)
#' @param prob.th numeric; probability threshold (default: 0.8)
#' @param TH.gene.exp.freq double; for a panel of signature genes, it will be classified as expressed if more than this value of genes are expressed (default: 0.65)
#' @param TH.gene.exp.freq.train.pos double; for a panel of signature genes, it will be classified as positive instance if more than this value of genes are expressed (default: 0.8)
#' @param TH.gene.exp.freq.train.neg double; for a panel of signature genes, it will be classified as negative instance if less than this value of genes are expressed (default: 0.2)
#' @param TH.silhouette double; threshold (default: -Inf)
#' @param RF.selectVar logical; (default: FALSE)
#' @param ... parameters passed to some method
#' @importFrom data.table as.data.table data.table `:=` melt dcast
#' @importFrom SummarizedExperiment assay rowData `colData<-`
#' @importFrom ggpubr ggscatter
#' @importFrom ggplot2 ggsave
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom plyr l_ply ldply
#' @importFrom S4Vectors DataFrame
#' @importFrom utils head
#' @details classify cells using the signature genes
#' @export
classifyCell.by.sigGene <- function(obj.list,gene.core.tb,pca.rotation=NULL,assay.name="exprs",out.prefix=NULL,
                                    adjB=NULL,meta.cluster=NULL,method="posFreq",
									gene.pos.tb=NULL,gene.neg.tb=NULL,
									meta.info.tb=NULL,meta.train.tb=NULL,ndownsample=1500,myseed=123456,
									verbose=F,
                                    sig.prevelance=0.65,ntop=10,ncores=16,
									bin.z.pos.th=0.3,bin.z.neg.th=0.3,
									prob.th=0.8,
									TH.gene.exp.freq=0.65,
									TH.gene.exp.freq.train.pos=0.8,TH.gene.exp.freq.train.neg=0.2,
									TH.silhouette=-Inf,
									RF.selectVar=F,...)
{
	if(is.null(meta.cluster)){
		mcls <- sort(unique(gene.core.tb$meta.cluster))
	}else{
		mcls <- sort(meta.cluster)
	}

	RhpcBLASctl::omp_set_num_threads(1)
	registerDoParallel(cores = ncores)

	getExpDataFromGeneSymbol <- function(obj,assay.name,gene.used,adjB=NULL){
		dat.block <- assay(obj,assay.name)[gene.used,,drop=F]
		#### adjB
		if(!is.null(adjB)){
		  dat.block <- simple.removeBatchEffect(dat.block,batch=obj[[adjB]])
		}
		dat.block <- t(scale(t(dat.block)))
		rownames(dat.block) <- rowData(obj)[rownames(dat.block),"display.name"]
		return(dat.block)
	}

	if(method=="posFreq"){
		binExp.tb <- as.data.table(ldply(seq_along(obj.list),function(i){
				  obj <- obj.list[[i]]
				  dataset.id <- names(obj.list)[i]
				  gene.used <- rownames(obj)[match(gene.core.tb$geneSymbol,rowData(obj)$display.name)]
				  gene.used <- unique(gene.used[!is.na(gene.used)])
				  dat.block <- getExpDataFromGeneSymbol(obj,assay.name,gene.used,adjB=adjB)
				  #### binarize all sig genes
				  dat.block.th <- (dat.block > bin.z.pos.th)
				  dat.block.bin <- dat.block.th

				  gene.pos.info <- NULL
				  gene.neg.info <- NULL
				  if(!is.null(gene.pos.tb)){
					  gene.pos.info <- llply(sort(unique(gene.pos.tb$meta.cluster)),function(xx){
								gene.xx <- gene.pos.tb[meta.cluster==xx,][["geneSymbol"]]
								gene.xx <- rownames(obj)[match(gene.xx,rowData(obj)$display.name)]
								gene.xx <- unique(gene.xx[!is.na(gene.xx)])
								dat.xx <- getExpDataFromGeneSymbol(obj,assay.name,gene.xx,adjB=adjB)
								dat.xx.th <- colSums(dat.xx > bin.z.pos.th)==length(gene.xx)
								names(dat.xx.th) <- colnames(dat.xx)
								return(dat.xx.th)
														})
					  names(gene.pos.info) <- sort(unique(gene.pos.tb$meta.cluster))
				  }
				  if(!is.null(gene.neg.tb)){
					  gene.neg.info <- llply(sort(unique(gene.neg.tb$meta.cluster)),function(xx){
								gene.xx <- gene.neg.tb[meta.cluster==xx,][["geneSymbol"]]
								gene.xx <- rownames(obj)[match(gene.xx,rowData(obj)$display.name)]
								gene.xx <- unique(gene.xx[!is.na(gene.xx)])
								dat.xx <- getExpDataFromGeneSymbol(obj,assay.name,gene.xx,adjB=adjB)
								dat.xx.th <- colSums(dat.xx < bin.z.neg.th)==length(gene.xx)
								names(dat.xx.th) <- colnames(dat.xx)
								return(dat.xx.th)
														})
					  names(gene.neg.info) <- sort(unique(gene.neg.tb$meta.cluster))
				  }

				  ### classification criteria
				  mcls.bin <- sapply(mcls,function(x){
							gene.sig <- gene.core.tb[meta.cluster==x,][["geneSymbol"]]
							exp.freq <- colSums(dat.block.bin[gene.sig,,drop=F])/length(gene.sig)
							ret <- exp.freq
							#ret <- structure(as.integer(exp.freq >= TH.gene.exp.freq),
							#				 names=colnames(dat.block.bin))
							return(ret)
							   })

				  binExp.tb.i <- cbind(data.table(cellID=rownames(mcls.bin),
												  dataset.id=names(obj.list)[i]),
									   mcls.bin)
				  idxCol.sigScore <- colnames(binExp.tb.i)[-c(1,2)]
				  binExp.tb.i$meta.cluster.top <- idxCol.sigScore[apply(binExp.tb.i[,idxCol.sigScore,with=F],1,
                                                            which.max)]
				  colnames(binExp.tb.i)[seq_along(idxCol.sigScore)+2] <- sprintf("score.%s",idxCol.sigScore)
				  for(xx in idxCol.sigScore){
					  binExp.tb.i[[sprintf("bin.%s",xx)]] <- as.integer(binExp.tb.i[[sprintf("score.%s",xx)]] >= TH.gene.exp.freq)
				  }
				  #### chose most likely and least likely cells for training
				  rf.tb <- as.data.table(ldply(idxCol.sigScore,function(xx){
					  cell.pos <- binExp.tb.i[[sprintf("score.%s",xx)]] >= TH.gene.exp.freq.train.pos
					  if(!is.null(gene.pos.info) && !is.null(gene.pos.info[[xx]])){
						  cell.pos <- cell.pos & gene.pos.info[[xx]]
					  }
					  if(!is.null(gene.neg.info) && !is.null(gene.neg.info[[xx]])){
						  cell.pos <- cell.pos & gene.neg.info[[xx]]
					  }
					  cell.pos <- as.integer(cell.pos)
					  cell.neg <- as.integer(binExp.tb.i[[sprintf("score.%s",xx)]] < TH.gene.exp.freq.train.neg)

					  if(sum(cell.neg)>sum(cell.pos)){
						  cell.neg <- sample(binExp.tb.i$cellID[which(cell.neg>0)],sum(cell.pos))
						  cell.pos <- binExp.tb.i$cellID[which(cell.pos>0)]
					  }else{
						  cell.pos <- sample(binExp.tb.i$cellID[which(cell.pos>0)],sum(cell.neg))
						  cell.neg <- binExp.tb.i$cellID[which(cell.neg>0)]
					  }
					  gene.train <- gene.core.tb[meta.cluster==xx,][["geneSymbol"]]
					  dat.train <- rbind(t(dat.block[gene.train,cell.pos,drop=F]),
										 t(dat.block[gene.train,cell.neg,drop=F]))
					  dat.train.label <- factor(c(rep("pos",length(cell.pos)),rep("neg",length(cell.neg))))
					  dat.pred <- t(dat.block[gene.train,,drop=F])
					  if(nrow(dat.train)>=20){
						  res.RF <-	run.RF(dat.train, dat.train.label, dat.pred, do.norm=F,
														 ntree = 500, ntreeIterat = 500,selectVar=RF.selectVar)
						  out.tb <- data.table(cellID=rownames(res.RF$yres),
											   meta.cluster=xx,
											   prob=res.RF$yres[,"pos"],
											   label=as.integer(res.RF$yres[,"pos"]>prob.th))
					  }else{
						  warning(sprintf("The number of instances for training is low (%s)\n",dataset.id))
						  out.tb <- data.table(cellID=rownames(dat.pred),
											   meta.cluster=xx,
											   prob=-1,
											   label=-1)
					  }
					  return(out.tb)
				  },.parallel=T))
				  rf.tb.prob <- dcast(rf.tb[,c("cellID","meta.cluster","prob")],cellID~meta.cluster,value.var="prob")
				  rf.tb.label <- dcast(rf.tb[,c("cellID","meta.cluster","label")],cellID~meta.cluster,value.var="label")
				  rf.tb.prob$RF.prob.max <- apply(rf.tb.prob[,idxCol.sigScore,with=F],1,
												  function(x){ x <- as.numeric(x); x[order(-x)[1]] })
				  rf.tb.prob$RF.prob.sec <- apply(rf.tb.prob[,idxCol.sigScore,with=F],1,
												  function(x){ x <- as.numeric(x); x[order(-x)[2]] })
				  rf.tb.prob$RF.meta.cluster.max <- apply(rf.tb.prob[,idxCol.sigScore,with=F],1,
												  function(x){ x <- as.numeric(x); idxCol.sigScore[order(-x)[1]] })
				  rf.tb.prob$RF.meta.cluster.sec <- apply(rf.tb.prob[,idxCol.sigScore,with=F],1,
												  function(x){ x <- as.numeric(x); idxCol.sigScore[order(-x)[2]] })
				  rf.tb.prob[,RF.class.max:=RF.meta.cluster.max]
				  rf.tb.prob[RF.prob.max<prob.th,RF.class.max:="unkown"]
				  rf.tb.prob[,RF.class.sec:=RF.meta.cluster.sec]
				  rf.tb.prob[RF.prob.sec<prob.th,RF.class.sec:="unkown"]
				  colnames(rf.tb.prob)[1+seq_along(idxCol.sigScore)] <- sprintf("RF.prob.%s",colnames(rf.tb.prob)[1+seq_along(idxCol.sigScore)])
				  colnames(rf.tb.label)[1+seq_along(idxCol.sigScore)] <- sprintf("RF.bin.%s",colnames(rf.tb.label)[1+seq_along(idxCol.sigScore)])
				  rf.tb <- merge(rf.tb.prob,rf.tb.label,by="cellID")
				  binExp.tb.i <- merge(binExp.tb.i,rf.tb,by="cellID")

				  if(verbose && !is.null(out.prefix)){
					l_ply(idxCol.sigScore,function(xx){
					  gene.plot <- gene.core.tb[meta.cluster==xx,][["geneID"]]
					  obj.tmp <- ssc.build(dat.block[gene.plot,binExp.tb.i$cellID,drop=F],assay.name="scaled.exprs")
					  cat(sprintf("all(colnames(obj.tmp)==binExp.tb.i$cellID):\n"))
					  print(all(colnames(obj.tmp)==binExp.tb.i$cellID))
					  colData(obj.tmp) <- DataFrame(binExp.tb.i)
					  colnames(obj.tmp) <- binExp.tb.i$cellID
					  p <- ssc.plot.violin(obj.tmp,"scaled.exprs",gene=gene.plot,
										   group.var=sprintf("RF.class.max"),clamp=c(0,6))
					  ggsave(sprintf("%s.RF.%s.%s.pdf",out.prefix,dataset.id,xx),p,width=7,height=9)
												  },.parallel=T)
				  }
				  return(binExp.tb.i)
							   }))

	}else if(method=="corrCtrl"){

		dat.plot.tb <- as.data.table(ldply(seq_along(mcls),function(j){
			gene.core.tb.j <- gene.core.tb[meta.cluster==mcls[j],]
			dat.plot.j <- ldply(seq_along(obj.list),function(i){
								  obj <- obj.list[[i]]
								  features <- list(gene.core.tb.j$geneID)
								  names(features) <- sprintf("sigScore.%s",mcls[j])
								  obj <- ssc.moduleScore(obj,features,assay.name="norm_exprs",adjB="batchV")
								  sig.score <- obj[[sprintf("sigScore.%s",mcls[j])]]
								  data.table(cellID=colnames(obj),
											 dataset.id=names(obj.list)[i],
											 meta.cluster=mcls[j],
											 sig.score=sig.score)
										},.parallel=T)
			return(dat.plot.j)
		}))

		binExp.tb <- dcast(dat.plot.tb,cellID+dataset.id~meta.cluster,value.var="sig.score")
		idxCol.sigScore <- colnames(binExp.tb)[-c(1,2)]
		binExp.tb$meta.cluster.top <- idxCol.sigScore[apply(binExp.tb[,idxCol.sigScore,with=F],1,
                                                            which.max)]
		colnames(binExp.tb)[seq_along(idxCol.sigScore)+2] <- sprintf("score.%s",idxCol.sigScore)
		for(xx in idxCol.sigScore){
			binExp.tb[[sprintf("bin.%s",xx)]] <- as.integer(binExp.tb[[sprintf("score.%s",xx)]] >= 1)
		}
	}else if(method=="pool" && !is.null(meta.info.tb)){

		exp.z.tb <- as.data.table(ldply(seq_along(obj.list),function(i){
				  obj <- obj.list[[i]]
				  dataset.id <- names(obj.list)[i]
				  if(is.null(pca.rotation)){
					  gene.used <- rownames(obj)[match(gene.core.tb$geneSymbol,rowData(obj)$display.name)]
				  }else{
					  gene.used <- rownames(obj)[match(rownames(pca.rotation),rowData(obj)$display.name)]
				  }
				  gene.used <- unique(gene.used[!is.na(gene.used)])
				  ### scale in the whole range of the dataset
				  dat.block <- getExpDataFromGeneSymbol(obj,assay.name,gene.used,adjB=adjB)
				  obj$ClusterID <- meta.info.tb$ClusterID[match(colnames(obj),meta.info.tb$cellID)]
				  obj$meta.cluster <- meta.info.tb$metaCluster[match(colnames(obj),meta.info.tb$cellID)]
				  obj <- ssc.downsample(obj,group.var="ClusterID",priority="silhouette",rd="seurat.pca")
				  ### select cells usd for training
				  #f.cell <- intersect(colnames(obj),meta.train.tb$cellID)
				  #obj <- obj[,f.cell]
				  #dat.block <- dat.block[,f.cell]
				  print(all(colnames(dat.block)==colnames(obj)))
				  ret.tb <- data.table(cellID=colnames(dat.block),
									   dataset.id=dataset.id,
									   ClusterID=obj$ClusterID,
									   meta.cluster=obj$meta.cluster,
									   silhouette=obj$silhouette,
									   silhouette.rank=obj$silhouette.rank)
				  if(is.null(pca.rotation)){
					  ret.tb <- cbind(ret.tb,t(dat.block))
				  }else{
					  dat.pc <- t(dat.block) %*% pca.rotation
					  ret.tb <- cbind(ret.tb,dat.pc)
				  }
				  ret.tb$cellID <- sprintf("%s.%s",dataset.id,ret.tb$cellID)
				  return(ret.tb)
		}))
		set.seed(myseed)
		exp.z.tb <- exp.z.tb[sample(nrow(exp.z.tb),replace=F),]
		dat.train.tb <- exp.z.tb[meta.cluster %in% unique(sort(gene.core.tb$meta.cluster)) &
									ClusterID %in% meta.train.tb$ClusterID,
								 ][silhouette>TH.silhouette,
								 ][,head(.SD,n=ndownsample),by="meta.cluster"]
		dat.test.tb <- exp.z.tb[meta.cluster %in% unique(sort(gene.core.tb$meta.cluster)) &
									ClusterID %in% meta.train.tb$ClusterID,
								 ][silhouette>TH.silhouette,
								 ][!(cellID %in% dat.train.tb$cellID),
								 ][,head(.SD,n=ndownsample/2),by="meta.cluster"]
		print(dat.train.tb[,.N,by=c("dataset.id","meta.cluster")])
		##exp.z.tb[1:4,1:7]
		##dat.train.tb[1:4,1:7]
		res.RF <-	run.RF(xdata=as.matrix(dat.train.tb[,-c(1:6)]),
									  xlabel=factor(dat.train.tb$meta.cluster),
									  ydata=as.matrix(exp.z.tb[,-c(1:6)]),
									  xdata.test=as.matrix(dat.test.tb[,-c(1:6)]),
									  xlabel.test=factor(dat.test.tb$meta.cluster),
									  do.norm=F,
									  selectVar=RF.selectVar,ncores=ncores,...)
		if(verbose){
			saveRDS(res.RF,file=sprintf("%s.res.RF.rds",out.prefix))
		}

		out.tb <- exp.z.tb[,c(1,4)]
		colnames(out.tb)[2] <- "meta.cluster.raw"
		out.tb$meta.cluster.RF <- res.RF$ylabel
		cls.set <- colnames(res.RF$yres)
		out.tb$meta.cluster.RF <- apply(res.RF$yres,1,function(x){ cls.set[which.max(x)] })
		out.tb$meta.cluster.RF.prob <- apply(res.RF$yres,1,function(x){ max(x) })
		out.tb$meta.cluster.RF.sec <- apply(res.RF$yres,1,function(x){ cls.set[order(-x)[2]] })
		out.tb$meta.cluster.RF.sec.prob <- apply(res.RF$yres,1, function(x){ x[order(-x)[2]] })
		prob.mtx <- res.RF$yres
		colnames(prob.mtx) <- sprintf("RF.prob.%s",colnames(prob.mtx))
		out.tb <- cbind(out.tb,prob.mtx)

		out.tb[,table(meta.cluster.RF,meta.cluster.raw)]
		meta.info.tb.a <- as.data.table(meta.info.tb)[,-c("metaCluster","metaCluster.old")]
		meta.info.tb.a[,cellID.raw:=sprintf("%s",cellID)]
		meta.info.tb.a[,cellID:=sprintf("%s.%s",dataset,cellID)]
		binExp.tb <- merge(meta.info.tb.a,out.tb,by="cellID")
		binExp.tb[,table(meta.cluster.raw,meta.cluster.RF)]

	}
	return(binExp.tb)
}


#' classify cells using pre-trained model
#' @param obj.list object; named list of object of \code{singleCellExperiment} class
#' @param mod.rf object; pre-trained model used for prediction
#' @param pca.rotation matrix; rotation matrix (default: NULL)
#' @param assay.name character; which assay (default: "exprs")
#' @param out.prefix character; output prefix (default: NULL).
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param meta.info.tb data.frame; cell information table (default: NULL)
#' @param pad.missing logical; (default: FALSE)
#' @param allZeroAsLow logical; (default: FALSE)
#' @param verbose logical; (default: FALSE)
#' @param ncores integer; number of CPU to used (default: 16)
#' @param ... parameters passed to some method
#' @importFrom data.table data.table `:=`
#' @importFrom  ggpubr ggscatter
#' @importFrom  SummarizedExperiment assay
#' @details classify cells using pre-trained model, which can be obtained by run classifyCell.by.sigGene()
#' @export
classifyCell.by.model <- function(obj.list,mod.rf,assay.name="exprs",out.prefix=NULL,
                                    adjB=NULL, meta.info.tb=NULL,pca.rotation=NULL,
									pad.missing=F,allZeroAsLow=F,
									verbose=F,ncores=16, ...)
{

	RhpcBLASctl::omp_set_num_threads(1)
	registerDoParallel(cores = ncores)

	getExpDataFromGeneSymbol <- function(obj,assay.name,gene.used,adjB=NULL,gene.pad=NULL,allZeroAsLow=F){
		gene.used <- intersect(gene.used,rownames(obj))
		dat.block <- assay(obj,assay.name)[gene.used,,drop=F]
		rownames(dat.block) <- rowData(obj)[rownames(dat.block),"display.name"]
		if(!is.null(gene.pad)){
			mtx.pad <- matrix(rep(0,length(gene.pad)*ncol(dat.block)),nrow=length(gene.pad))
			rownames(mtx.pad) <- gene.pad
			dat.block <- rbind(dat.block,mtx.pad)
		}
		#### adjB
		if(!is.null(adjB)){
			dat.block <- simple.removeBatchEffect(dat.block,batch=obj[[adjB]])
		}
		if(allZeroAsLow){
			f.gene <- rowSums(dat.block==0)==ncol(dat.block)
			dat.block.a <- t(scale(t(dat.block[!f.gene,,drop=F])))
			#### z-score 10000 is very low
			if(sum(f.gene)>0)
			{
				dat.block.b <- matrix(rep(-10000,sum(f.gene)*ncol(dat.block)),nrow=sum(f.gene))
				rownames(dat.block.b) <- rownames(dat.block)[f.gene]
				dat.block <- rbind(dat.block.a,dat.block.b)
			}else{
				dat.block <- dat.block.a
			}
		}else{
			dat.block <- t(scale(t(dat.block)))
		}
		return(dat.block)
	}
	

	{

		exp.z.tb <- as.data.table(ldply(seq_along(obj.list),function(i){
				  obj <- obj.list[[i]]
				  dataset.id <- names(obj.list)[i]
				  gene.pad <- NULL
				  if(is.null(pca.rotation)){
					  ###gene.used <- rownames(obj)[match(gene.core.tb$geneSymbol,rowData(obj)$display.name)]
					  gene.used <- rownames(obj)[match(rownames(mod.rf$importance),rowData(obj)$display.name)]
					  if(any(is.na(gene.used))){
						  if(!pad.missing){
							  cat(sprintf("Not all gene(s) found in dataset %s \n",dataset.id))
							  return(NULL)
						  }else{
							  idx.na <- which(is.na(gene.used))
							  gene.pad <- rownames(mod.rf$importance)[idx.na]
						  }
					  }
				  }else{
					  gene.used <- rownames(obj)[match(rownames(pca.rotation),rowData(obj)$display.name)]
				  }
				  
				  ###cat(sprintf("dataset: %s\n",dataset.id))
				  dat.block <- getExpDataFromGeneSymbol(obj,assay.name,gene.used,adjB=adjB,
														gene.pad=gene.pad,allZeroAsLow=allZeroAsLow)
				  if(sum(is.na(dat.block)) > 0){
					  cat(sprintf("NA value(s) found in dataset %s (may casuesed by all 0s in some genes)\n",dataset.id))
					  return(NULL)
				  }
				  dat.block <- dat.block[rownames(mod.rf$importance),]

				  ret.tb <- data.table(cellID=colnames(dat.block))
				  #ret.tb$meta.cluster=meta.info.tb$metaCluster[match(ret.tb$cellID,meta.info.tb$cellID)]
				  #ret.tb$ClusterID=meta.info.tb$ClusterID[match(ret.tb$cellID,meta.info.tb$cellID)]
				  if(is.null(pca.rotation)){
					  ret.tb <- cbind(ret.tb,t(dat.block))
				  }else{
					  dat.pc <- t(dat.block) %*% pca.rotation
					  ret.tb <- cbind(ret.tb,dat.pc)
				  }
				  ret.tb$cellID <- sprintf("%s.%s",dataset.id,ret.tb$cellID)
				  return(ret.tb)
		}))

		##### predict
		ydata <- as.matrix(exp.z.tb[,-c(1)])
		rownames(ydata) <- exp.z.tb$cellID
		#### ranger:::predict.ranger
		if("ranger" %in% class(mod.rf)){
			yres <- predict(mod.rf,ydata,num.threads=ncores)$predictions
		}else{
			yres <- predict(mod.rf,ydata,type="prob")
		}
		cls.set <- colnames(yres)
		ylabel <- apply(yres,1,function(x){ cls.set[which.max(x)] })
		names(ylabel) <- rownames(ydata)

	    ### output table
		out.tb <- data.table(cellID=exp.z.tb$cellID)
		out.tb$meta.cluster.RF <- ylabel
		out.tb$meta.cluster.RF.prob <- apply(yres,1,function(x){ max(x) })
		out.tb$meta.cluster.RF.sec <- apply(yres,1,function(x){ cls.set[order(-x)[2]] })
		out.tb$meta.cluster.RF.sec.prob <- apply(yres,1, function(x){ x[order(-x)[2]] })
		prob.mtx <- yres
		colnames(prob.mtx) <- sprintf("RF.prob.%s",colnames(prob.mtx))
		out.tb <- cbind(out.tb,prob.mtx)

		if(!is.null(meta.info.tb)){
			meta.info.tb.a <- as.data.table(meta.info.tb)
			meta.info.tb.a[,cellID.raw:=sprintf("%s",cellID)]
			meta.info.tb.a[,cellID:=sprintf("%s.%s",dataset,cellID)]
			meta.info.tb.a[,meta.cluster.raw:=""]
			meta.info.tb.b <- merge(meta.info.tb.a,out.tb,by="cellID")
		}else{
			meta.info.tb.b <- out.tb
		}

	}
	return(meta.info.tb.b)
}


#' plot signature genes
#' @param obj.list object; named list of object of \code{singleCellExperiment} class
#' @param gene.core.tb data.frame; signature genes 
#' @param out.prefix character; output prefix 
#' @param meta.info.tb data.frame; cell information table (default: NULL)
#' @param meta.col character; column in meta.info.tb, by which to group cells (default: "meta.cluster")
#' @param assay.name character; which assay (default: "exprs")
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param meta.cluster character; (default: NULL)
#' @param verbose logical; (default: FALSE)
#' @param ptype character; (default: "heatmap")
#' @param ... parameters passed to plotting methods
#' @param ncores integer; number of CPU to used (default: 16)
#' @param p.ncol integer; number of columns in the violin plot (default: 2)
#' @param pdf.width double; width of the output plot (default: NULL)
#' @param pdf.height double; height of the output plot (default: NULL)
#' @importFrom data.table data.table melt dcast
#' @importFrom ggplot2 ggsave
#' @importFrom ggpubr ggscatter
#' @importFrom grid unit gpar
#' @importFrom SummarizedExperiment assay
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @details make some plot(s) of signature genes
#' @export
plotSigGene <- function(obj.list,gene.core.tb,out.prefix,assay.name="exprs",
                                    adjB=NULL,meta.cluster=NULL,
									meta.info.tb=NULL,meta.col="metaCluster",
									p.ncol=2,pdf.width=NULL,pdf.height=NULL,
									verbose=F,ptype="heatmap", ncores=16,...)
{
	if(is.null(meta.cluster)){
		mcls <- sort(unique(gene.core.tb$meta.cluster))
	}else{
		mcls <- sort(meta.cluster)
	}

	RhpcBLASctl::omp_set_num_threads(1)
	registerDoParallel(cores = ncores)

	getExpDataFromGeneSymbol <- function(obj,assay.name,gene.used,adjB=NULL){
		dat.block <- assay(obj,assay.name)[gene.used,,drop=F]
		#### adjB
		if(!is.null(adjB)){
		  dat.block <- simple.removeBatchEffect(dat.block,batch=obj[[adjB]])
		}
		dat.block <- t(scale(t(dat.block)))
		rownames(dat.block) <- rowData(obj)[rownames(dat.block),"display.name"]
		return(dat.block)
	}
	
	if(ptype=="heatmap"){
		exp.z.tb <- as.data.table(ldply(seq_along(obj.list),function(i){
				  obj <- obj.list[[i]]
				  dataset.id <- names(obj.list)[i]
				  gene.used <- rownames(obj)[match(gene.core.tb$geneSymbol,rowData(obj)$display.name)]
				  gene.used <- unique(gene.used[!is.na(gene.used)])
				  dat.block <- getExpDataFromGeneSymbol(obj,assay.name,gene.used,adjB=adjB)
				  ret.tb <- data.table(cellID=colnames(dat.block))
				  ret.tb$meta.cluster=meta.info.tb[[meta.col]][match(ret.tb$cellID,meta.info.tb$cellID)]
				  ret.tb$ClusterID=meta.info.tb$ClusterID[match(ret.tb$cellID,meta.info.tb$cellID)]
				  ret.tb <- cbind(ret.tb,t(dat.block))
				  ret.tb$cellID <- sprintf("%s.%s",dataset.id,ret.tb$cellID)
				  return(ret.tb)
		}))
		###exp.z.tb.debug <<- exp.z.tb
		exp.z.ave.tb <- melt(exp.z.tb[,-c("cellID")],
				   id.vars=c("meta.cluster","ClusterID"))[,.(value=mean(value)),by=c("meta.cluster","ClusterID","variable")]
		exp.z.ave.tb <- dcast(exp.z.ave.tb,meta.cluster+ClusterID~variable)
		###exp.z.ave.tb[1:4,1:5]

		sce.tmp <- ssc.build(t(exp.z.ave.tb[,-c("meta.cluster","ClusterID")]))
		colnames(sce.tmp) <- exp.z.ave.tb$ClusterID
		sce.tmp$ClusterID <- exp.z.ave.tb$ClusterID
		##sce.tmp$meta.cluster <- exp.z.ave.tb$meta.cluster
		##colData(sce.tmp)[[meta.col]] <- exp.z.ave.tb$meta.cluster
		sce.tmp[[meta.col]] <- exp.z.ave.tb$meta.cluster
		gene.core.plot.tb <- gene.core.tb[!duplicated(geneID),]
		o.width <- 20
		o.height <- 10
		if(!is.null(pdf.width)){ o.width <- pdf.width }
		if(!is.null(pdf.height)){ o.height <- pdf.height }
		ssc.plot.heatmap(sce.tmp[,!grepl("^unk",sce.tmp[[meta.col]])],
				 out.prefix=sprintf("%s.gene.core.zscore.%s",out.prefix,meta.col),
				 columns=meta.col,columns.order=meta.col,
				 gene.desc=gene.core.plot.tb,
				 row.split=gene.core.plot.tb$Group,
				 column.split=sce.tmp[[meta.col]],
				 row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
				 row_title_rot = 0, column_title_gp = gpar(fontsize = 0),
				 border = TRUE,
				 pdf.width=o.width,pdf.height=o.height,do.scale=F,
				 z.lo=-0.8,z.hi=0.8,z.step=0.2,exp.title="Z-Score",
				 do.clustering.row=F,
				 do.clustering.col=T,...)
	}else if(ptype=="violin"){
		o.width <- 30
		o.height <- 9
		if(!is.null(pdf.width)){ o.width <- pdf.width }
		if(!is.null(pdf.height)){ o.height <- pdf.height }
		l_ply(seq_along(obj.list),function(i){
			obj <- obj.list[[i]]
			dataset.id <- names(obj.list)[i]
			loginfo(sprintf("plotting violin (%s):",dataset.id))
			l_ply(mcls,function(xx){
				  gene.plot <- gene.core.tb[meta.cluster==xx,][["geneSymbol"]]
				  gene.used <- rownames(obj)[match(gene.plot,rowData(obj)$display.name)]
				  gene.used <- unique(gene.used[!is.na(gene.used)])
				  dat.block <- getExpDataFromGeneSymbol(obj,assay.name,gene.used,adjB=adjB)
				  ##colnames(dat.block) <- sprintf("%s.%s",dataset.id,colnames(dat.block))
				  obj.tmp <- ssc.build(dat.block[,,drop=F],assay.name="scaled.exprs")
				  obj.tmp[[meta.col]] <- meta.info.tb[[meta.col]][match(colnames(obj.tmp),
																		meta.info.tb$cellID)]
				  p <- ssc.plot.violin(obj.tmp,"scaled.exprs",gene=rownames(obj.tmp),
									   group.var=sprintf(meta.col),clamp=c(0,6),p.ncol=p.ncol,...)
				  ggsave(sprintf("%s.RF.%s.%s.pdf",out.prefix,dataset.id,xx),p,width=o.width,height=o.height)
											  },.parallel=T)
								  })
	}
}

#' binarize the expression value using the distribution
#' @param x object; vector
#' @param out.prefix character; output prefix [default: NULL]
#' @param G integer; number of components [default: NULL]
#' @param topNAsHi integer; treat this number of top components as high;if it's zero, components with mean larger than 0 are high [default: 1]
#' @param e.TH double; for plot purpose: comparing the components with the user defined threhold [default: NULL]
#' @param e.name character; name of the expression metric [default: "Exp"]
#' @param verbose logical; [default: F]
#' @param draw.CI logical; [default: T]
#' @param zero.as.low logical; [default: T]
#' @param run.extremevalue logical; [default: F]
#' @param zero.notUsed logical; [default: F]
#' @param my.seed integer; [default: 9997]
#' @param ... parameters passed to plot.densityMclust
#' @importFrom data.table data.table
#' @importFrom mclust densityMclust
#' @importFrom graphics abline layout plot lines rect
#' @importFrom grDevices pdf png dev.off
#' @importFrom stats ppoints lm qnorm dnorm
#' @details use mixture model to classify each data point. If G is null, the largest component is corresponding to binarized expression 1, and other components are corresponding to binarized expressed 0.
#' @export
binarizeExp <- function(x,out.prefix=NULL,G=NULL,topNAsHi=1,e.TH=NULL,e.name="Exp",verbose=F,
                         draw.CI=T, zero.as.low=T,run.extremevalue=F,zero.notUsed=F,my.seed=9997,...)
{
  ###require(mclust)
  if(!is.null(names(x))){
    o.df <- data.frame(sample=names(x),bin.Exp=-1,o.Exp=x,stringsAsFactors = F)
  }else{
    names(x) <- sprintf("S%04d",seq_along(x))
    o.df <- data.frame(sample=sprintf("S%04d",seq_along(x)),bin.Exp=-1,o.Exp=x,stringsAsFactors = F)
  }
  rownames(o.df) <- o.df$sample
  if(zero.notUsed){
    f<-is.finite(x) & x>0
  }else{
	f<-is.finite(x)
  }
  if(sum(f)<3){
      colnames(o.df) <- c("sample",e.name)
      if(verbose){
        return(list(x_mix=NULL,o.df=o.df))
      }else{
        return(structure(o.df[,e.name],names=rownames(o.df)))
      }
  }

  x <- x[f]
  set.seed(my.seed)
  x_mix <- mclust::densityMclust(x,G=G,modelNames=c("E","V"))
  x_mix_summary <- summary(x_mix)

  quantEstDist <- mclust::quantileMclust(x_mix, p = ppoints(length(x)))
  quantSample <- sort(x_mix$data)
  #quantEstDist.debug <<- quantEstDist
  #x.debug <<- x
  #x_mix.debug <<- x_mix
  ret.gof <- list(R2=NA)
  if(verbose){
	  tryCatch({
		  res.lm <- lm(y~x,data=data.frame(x=quantEstDist,y=quantSample))
		  res.lm.summ <- summary(res.lm)
		  ret.gof <- list(R2=res.lm.summ$r.squared)
	  },error=function(e){
		  print(e)
	  })
  }

  dat.extra <- NULL
  if(run.extremevalue){
    ###dat.extra <- classify.outlier(x,out.prefix=sprintf("%s.extremevalue",out.prefix))
    dat.extra <- classify.outlier(x,out.prefix=NULL)
  }

  if(verbose && !is.null(out.prefix)){
	  print(x_mix_summary)
      if(grepl(".png$",out.prefix,perl = T)){
          png(out.prefix,width=800,height=600)
	      old_par<-par(no.readonly=T)
	      layout(matrix(c(1,1,1,1,1,1,2,3,4), 3, 3, byrow = TRUE))
	      a_par<-par(cex.axis=2,cex.lab=2,cex.main=1.8,mar=c(5,6,4,2)+0.1)
      }else{
          pdf(out.prefix,width=8,height=6)
      }
	  plot(x_mix,what="density",data=x,breaks=50,col="darkgreen",lwd=2,main="",...)
	  abline(v=x_mix_summary$mean,lty=2)
	  if(!is.null(e.TH)){
		abline(v=e.TH,lty=2,col="red")
	  }

	  for(i in 1:x_mix_summary$G)
	  {
		i_mean<-x_mix_summary$mean[i]
		i_sd<-sqrt(x_mix_summary$variance[i])
		i_pro<-x_mix_summary$pro[i]
		#i_sd<-RC_mix_summary$variance[i]
		d<-qnorm(c(0.0013,0.9987),i_mean,i_sd)
		e<-i_pro*dnorm(i_mean,i_mean,i_sd)
		lines(seq(d[1],d[2],by=0.01),i_pro*dnorm(seq(d[1],d[2],by=0.01),i_mean,i_sd),col="orange",lwd=2)
        if(draw.CI){ rect(d[1],0,d[2],e+0.02,col=NA,border="blue") }
		#rect(d[1],0,d[2],e+0.02,col=rgb(0,0,0.8,0.2),border=NA)
	  }

	  plot(x_mix,data=x,breaks=20,col="darkgreen",lwd=2,what="BIC")
	  mclust::densityMclust.diagnostic(x_mix,type = "cdf",cex.lab=1.5)
	  tryCatch({
		mclust::densityMclust.diagnostic(x_mix,type = "qq")
	  },error=function(e){ cat(sprintf("Error in  densityMclust.diagnostic(x_mix,type = \"qq\")\n")); print(e); e })
	  dev.off()
  }

  o.df[names(x_mix$classification),"classification"] <- x_mix$classification
  o.df[names(x_mix$classification),"bin.Exp"] <- x_mix$classification
  colnames(o.df)[2] <- e.name
  if(is.null(G)){
	  if(topNAsHi==0){
		iG.2nd <- x_mix$G - sum(x_mix_summary$mean > 0)
	  }else{
		iG.2nd <- x_mix$G-topNAsHi
	  }
      f.low <- o.df[[e.name]] <= iG.2nd
      o.df[[e.name]][f.low] <- 0
      o.df[[e.name]][!f.low] <- 1
  }else if(!is.null(G) && x_mix_summary$G==3){
	  ### determin which classes 'not-expressed' and which classes 'expressed'
	  i_mean<-x_mix_summary$mean
	  i_sd<-sqrt(x_mix_summary$variance)
	  ci95.1 <- qnorm(c(0.0013,0.9987),i_mean[1],i_sd[1])
	  ci95.1.overlap <- pnorm(ci95.1[2],i_mean[2],i_sd[2])
	  if(ci95.1.overlap>0.3333){
		  C.min <-  2
	  }else{
		  C.min <- 1
	  }
	  ci95.2 <- qnorm(c(0.0013,0.9987),i_mean[2],i_sd[2])
	  ci95.3 <- qnorm(c(0.0013,0.9987),i_mean[3],i_sd[3])
	  ci95.3.overlap <- pnorm(ci95.3[1],i_mean[2],i_sd[2],lower.tail = F)
	  if(ci95.3.overlap>0.3333){
		  C.max <- 2
	  }else{
		  C.max <- 3
	  }
#	  ci95.2.overlap <- pnorm(ci95.2[2],i_mean[3],i_sd[3])
#	  if(ci95.2.overlap>0.3333){
#	  	  C.max <- 2
#	  }else{
#	  	  C.max <- 3
#	  }
      if(zero.as.low){
	    f.low <- o.df[,e.name] <= C.min
      }else{
	    f.low <- o.df[,e.name] <= C.min & o.df[,e.name] != -1
      }
	  f.hi <-  o.df[,e.name] >= C.max
	  f.mid <- (o.df[,e.name] > C.min) & (o.df[,e.name] < C.max)
	  o.df[f.low,e.name] <- 0
	  o.df[f.mid,e.name] <- 0.5
	  o.df[f.hi,e.name] <- 1
	  #f.G <- o.df[,e.name] < C.max
	  #o.df[f.G,e.name] <- 0
	  #o.df[!f.G,e.name] <- 1
  }else if(x_mix_summary$G==2){
      if(zero.as.low){
	    f.low <- o.df[,e.name] <= 1
      }else{
	    f.low <- o.df[,e.name] <= 1 & o.df[,e.name] != -1
      }
	  f.hi <-  o.df[,e.name] == 2
	  o.df[f.low,e.name] <- 0
	  o.df[f.hi,e.name] <- 1
  }
  if(verbose){
	return(list(x_mix=x_mix,o.df=o.df,gof=ret.gof,dat.extra=dat.extra))
  }else{
	return(structure(o.df[,e.name],names=rownames(o.df)))
  }
}


#' outlier detection using extremevalues
#' @param x object; vector
#' @param out.prefix character; output prefix [default: NULL]
#' @param e.name character; name of the expression metric [default: "Exp"]
#' @param th.score numeric; numeic vector to show the threshold specified manually [default: NULL]
#' @importFrom extremevalues outlierPlot getOutliers
#' @importFrom data.table as.data.table
#' @importFrom grDevices png pdf dev.off
#' @importFrom graphics par
#' @importFrom ggplot2 ggplot geom_line geom_vline ggsave aes aes_string geom_density theme_bw
#' @importFrom stats dnorm
#' @details outlier detection using extremevalues
#' @export
classify.outlier <- function(x,out.prefix=NULL,e.name="Exp",th.score=NULL)
{
	##### outlier detection (use method I at last)
	K <- extremevalues::getOutliers(x,method="I",distribution="normal")
	L <- extremevalues::getOutliers(x,method="II",distribution="normal")

	ii <- seq(K$limit[1], K$limit[2],0.01)
    x.fit <- dnorm(ii,mean = K$mu,sd = K$sigma)

    if(is.null(names(x))){ names(x) <- sprintf("S%04d",seq_along(x)) }
    ret.1 <- data.frame(sample=names(x),bin.Exp=0,o.Exp=x,stringsAsFactors=F)
    ##ret.1 <- data.table(score=x)
    ##ret.1[,score.cls:=0]
    ret.1$bin.Exp[K$iRight] <- 1
    ret.1$classification <- ret.1$bin.Exp
    colnames(ret.1)[2] <- e.name

    ret.2 <- as.data.table(list("score"=ii,"density"=x.fit))

	##### plot #####
    if(!is.null(out.prefix)){
        ##pdf(sprintf("%s.extremevalues.outlier.pdf",out.prefix),width = 10,height = 6)
        png(sprintf("%s.extremevalues.outlier.png",out.prefix),width = 1000,height = 600)
        opar <- par(mfrow=c(1,2))
        extremevalues::outlierPlot(x,K,mode="qq")
        extremevalues::outlierPlot(x,L,mode="residual")
        dev.off()
        par(opar)

	    p <- ggplot(ret.1, aes(o.Exp)) + geom_density(colour="black") + theme_bw() +
		        geom_line(data=ret.2,aes_string(x="score",y="density"),colour="red") +
		        geom_vline(xintercept = K$limit,linetype=2,colour="orange")
	    if(!is.null(th.score)){
			p <- p + geom_vline(xintercept = th.score,linetype=2,colour="black")
		}
	    ggsave(filename = sprintf("%s.extremevalues.density.pdf",out.prefix),width = 4,height = 3)
    }
    return(list("score.cls.tb"=ret.1,"fit.value.tb"=ret.2,"K"=K,"R2"=K$R2))
}


#' plot genes expression in pairs of clusters to examine the correlation
#' @param obj object; object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param x.cluster character; cluster to be showed in x axis
#' @param out.prefix character; output prefix.
#' @param genes vector; genes to used to calculate the correlation and make the scatter plot (default: NULL)
#' @param gene.highlight vector; genes to highlight (default: NULL)
#' @param plot.ncol integer; number of columns in facet_wrap() (default: 7)
#' @param plot.width double; width of the output plot (default: 22)
#' @param plot.height double; height of the output plot (default: 22)
#' @importFrom data.table setDT melt dcast
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 geom_hline geom_vline facet_wrap ggsave
#' @importFrom ggpubr ggscatter
#' @details used to examine the correlation between pairs of clusters
#' @export
plotPairsCor <- function(obj,assay.name,x.cluster,out.prefix,
                           genes=NULL,gene.highlight=NULL,
                           plot.ncol=7,plot.width=22,plot.height=22)
{
    #require("data.table")
    #require("ggplot2")
    #require("ggpubr")
    dat.assay <- assay(obj,assay.name)
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


