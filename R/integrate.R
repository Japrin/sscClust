#' integrate dataset by calculating the average expression and then clustering the pseudo bulk data
#' @param sce.list list; list of object of \code{singleCellExperiment} class
#' @param out.prefix character; output prefix.
#' @param assay.name character; which assay (default: "exprs")
#' @param is.avg logical; whether the objects in sce.list are average expression (default: FALSE)
#' @param ncores integer; number of cores to use (default: 6)
#' @param use.deg logical; whether use only the differentially expressed genes (default: FALSE)
#' @param gene.de.list list; if not NULL, each element is a data.frame with "geneID" column (default: NULL)
#' @param avg.by character; calculate the average expression of cells group by the specifid column (default: "majorCluster")
#' @param method.avg character; method of calculate the average expression. Passed to `avg` of `ssc.average.cell`.(default: "zscore")
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
    gene.common <- c()
    for(i in seq_along(sce.list)){
      if(length(gene.common)==0)
            gene.common <- rownames(sce.list[[i]])
        else
            gene.common <- intersect(gene.common,rownames(sce.list[[i]]))
    }

    ##Differential expressed genes
    if(use.deg){
        if(is.null(gene.de.list)){
            gene.de.list <- list()
            for(i in seq_along(sce.list)){
                de.out <- ssc.clusterMarkerGene(sce.list[[i]],n.cores=ncores,...)
                ###de.out <- ssc.clusterMarkerGene(sce.list[[i]],n.cores=ncores)
                gene.de.list[[i]] <- de.out$gene.table
            }
        }
        gene.de.common <- c()
        for(i in seq_along(gene.de.list)){
            if(length(gene.de.common)==0)
                gene.de.common <- gene.de.list[[i]]$geneID
            else
                gene.de.common <- c(gene.de.common,gene.de.list[[i]]$geneID)
        }
        gene.de.common <- unique(gene.de.common)
        gene.common <- intersect(gene.common,gene.de.common)
    }

    ### cal the average
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
                             return(obj.avg)
    },.parallel = F)
    names(sce.avg.list) <- names(sce.list)

    dat.avg.mtx <- NULL
    for(aid in names(sce.avg.list)){
        if(is.null(dat.avg.mtx)){
            dat.avg.mtx <- assay(sce.avg.list[[aid]],assay.name)
        }else{
            dat.avg.mtx <- cbind(dat.avg.mtx,assay(sce.avg.list[[aid]],assay.name))
        }
    }

    #### correlation between clusters
    dat.cor.mtx <- cor(dat.avg.mtx,method="pearson")

    ###  clustering the pseudo bulk data
    #all(rownames(dat.avg.mtx)==rownames(sce.avg.list[[1]]))
    sce.pb <- ssc.build(dat.avg.mtx,display.name=rowData(sce.avg.list[[1]])$display.name)

    rowData(sce.pb)$gene.common <- rownames(sce.pb) %in% gene.common
    sce.pb <- ssc.reduceDim(sce.pb,method="pca",method.vgene="gene.common",
                            pca.npc=15,seed=9997)

    #ssc.plot.pca(sce.pb)
    m <- regexec("^(.+?)\\.[A-Z0-9]",colnames(sce.pb),perl=T)
    mm <- regmatches(colnames(sce.pb),m)
    colData(sce.pb)[,avg.by] <- colnames(sce.pb)
    colData(sce.pb)[,"dataset.id"] <- sapply(mm,"[",2)
    colData(sce.pb)

    #### clustering the clusters
    #sce.pb <- ssc.clust(sce.pb, method.reduction="pca",
    #                method="hclust", k.batch=c(0),k.max=40,
    #                nboot=1000,seed = 9997)
    sce.pb <- ssc.clust(sce.pb, method.reduction="pca",
                    method="dynamicTreeCut",
                    deepSplit=4, minClusterSize=2, seed = 9997)

    p <- ssc.plot.tsne(sce.pb,columns = "dataset.id",reduced.name = "pca.tsne",size=3)
    ggsave(sprintf("%s.pca.tsne.aid.pdf",out.prefix),width=5,height=4)

    #p <- ssc.plot.tsne(sce.pb,columns = "pca.hclust.k0",
    #               reduced.name = "pca.tsne",
    #               size=3)
    #ggsave(sprintf("%s.pca.tsne.hclust.k0.pdf",out.prefix),width=5.5,height=4)

    p <- ssc.plot.tsne(sce.pb,columns = "pca.dynamicTreeCut.kauto",
                   reduced.name = "pca.tsne",
                   size=3)
    ggsave(sprintf("%s.pca.tsne.dynamicTreeCut.kauto.pdf",out.prefix),
           width=6.0,height=4)

    branch.out <- plot.branch(metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$hclust,
            cluster=sce.pb$pca.dynamicTreeCut.kauto,
            out.prefix=sprintf("%s.dynamicTreeCut",out.prefix))
    metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch <- branch.out$branch

    for(z.lo in c(0,0.5,-0.5,-1)){
        plot.matrix.simple(dat.cor.mtx,
                           out.prefix=sprintf("%s.avg.%s.cor.zlo.%s",out.prefix,"pearson",z.lo),
                           do.clust=metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch,
                           show.number = F,
                           z.lo = z.lo, z.hi = 1,exp.name="Cor")
    }

    plot.matrix.simple(as.matrix(metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$dist),
                       out.prefix=sprintf("%s.avg.pca.dist",out.prefix),
                       do.clust=metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch,show.number = F,
                       palatte=(brewer.pal(n = 7,name = "RdYlBu")),
                       z.lo = 0, z.hi = 15,exp.name="dist")

    #### heatmap show average expression of specified genes
    ###return(list("sce.pb"=sce.pb,"branch.out"=branch.out))
    if(use.deg && is.null(gene.de.list)){
        metadata(sce.pb)$ssc$gene.de.list <- gene.de.list
    }
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


