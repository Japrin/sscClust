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
							   z.lo=-15,z.hi=15,z.step=3,
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
	loginfo(sprintf("integrate.by.avg run successfully"))
    return(sce.pb)
}

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


