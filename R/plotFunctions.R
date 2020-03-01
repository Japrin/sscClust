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
#' @importFrom ggplot2 ggplot ggsave scale_colour_gradientn geom_point facet_wrap theme_bw coord_cartesian
#' @importFrom data.table melt data.table
#' @importFrom utils head
#' @param Y matrix or data.frame; Gene expression data, rownames shoud be gene id, colnames
#' should be sample id
#' @param dat.map data.frame; tSNE map, must be two columns data.frame and rownames should be sample id
#' @param gene.to.show character; gene id to be showd on the tSNE map
#' @param out.prefix character; output prefix (default: NULL)
#' @param p.ncol integer; number of columns in the plot's layout (default: 3)
#' @param theme.use function; which theme to use (default: theme_bw)
#' @param xlim integer or NULL; only draw points lie in the ragne specified by xlim and ylim (default NULL)
#' @param ylim integer or NULL; only draw points lie in the ragne specified by xlim and ylim (default NULL)
#' @param size double; points' size. If NULL, infer from number of points (default NULL)
#' @param width numeric; width of the plot (default: 9)
#' @param height numeric; height of the plot (default: 8)
#' @param pt.alpha numeric; alpha of the points (default: 0.5)
#' @param pt.order character; (default: "value")
#' @param clamp integer vector; expression values will be clamped to the range defined by this parameter, such as c(0,15). (default: NULL )
#' @param scales character; whether use the same scale across genes. one of "fixed" or "free" (default: "fixed")
#' @param vector.friendly logical; output vector friendly figure (default: FALSE)
#' @param par.legend list; lengend parameters, used to overwrite the default setting; (default: list())
#' @details For genes contained in both `Y` and `gene.to.show`, show their expression on the tSNE
#' map provided as `dat.map`. One point in the map represent a cell; cells with higher expression
#' also have darker color.
#' @return a ggplot object
ggGeneOnTSNE <- function(Y,dat.map,gene.to.show,out.prefix=NULL,p.ncol=3,theme.use=theme_bw,
                         xlim=NULL,ylim=NULL,size=NULL,pt.alpha=0.5,pt.order="value",clamp=NULL,
                         width=9,height=8,scales="fixed",vector.friendly=F,par.legend=list()){
  #suppressPackageStartupMessages(require("data.table"))
  #requireNamespace("ggplot2",quietly = T)
  #requireNamespace("RColorBrewer",quietly = T)
  if(!is.null(out.prefix)){
    dir.create(dirname(out.prefix),showWarnings = F,recursive = T)
  }
  f.g <- gene.to.show %in% rownames(Y)
  if(sum(!f.g)>0){
    warning(sprintf("Some genes not in the expression data: \n"))
    print(head(gene.to.show[!f.g]))
  }
  gene.to.show <- gene.to.show[f.g]

  dat.plot <- data.table::data.table(sample=rownames(dat.map),stringsAsFactors = F)
  dat.plot <- cbind(dat.plot,dat.map,t(as.matrix(Y[gene.to.show,dat.plot$sample,drop=F])))
  colnames(dat.plot) <- c("sample","Dim1","Dim2",names(gene.to.show))
  dat.plot.melt <- data.table::melt(dat.plot,id.vars = c("sample","Dim1","Dim2"))
  dat.plot.melt <- dat.plot.melt[!is.na(value),]
  if(pt.order=="value"){
	  dat.plot.melt <- dat.plot.melt[order(dat.plot.melt$value,decreasing = F),]
  }else if(pt.order=="random"){
	  dat.plot.melt <- dat.plot.melt[sample(nrow(dat.plot.melt),nrow(dat.plot.melt)),]
  }
  npts <- nrow(dat.plot)
  dat.plot.melt[,variable:=factor(variable,levels=names(gene.to.show),ordered=T)]
  
  make.plot <- function(gene.to.show,dat.plot.melt,show.legend=T,clamp="none",value.range=NULL,size=NULL,vector.friendly=F,out.prefix=NULL)
  {
	  lapply(names(gene.to.show),function(x){
						.dd <- dat.plot.melt[variable==x,]
						if(!is.null(clamp) && clamp=="none"){
						}else{
							if(is.null(clamp)){
								clamp <- quantile(.dd$value,c(0.05,0.95))
							}
							.dd[value < clamp[1],value:=clamp[1]]
							.dd[value > clamp[2],value:=clamp[2]]
						}
						if(is.null(value.range)){
							value.range <- pretty(.dd$value)
						}
						p <- ggplot2::ggplot(.dd,aes(Dim1,Dim2))+
							geom_point(aes(colour=value),
									   size=if(is.null(size)) sscClust:::auto.point.size(npts)*1.1 else size,
									   alpha=pt.alpha,stroke=0,shape=16) +
							labs(title=x, x ="", y = "")
						p <- p + do.call(scale_colour_gradientn,c(list(colours=RColorBrewer::brewer.pal(9,"YlOrRd"),
																	 limits=c(value.range[1],
																			  value.range[length(value.range)])),
																  par.legend))
						p <- p + theme.use()
					    p <- p + theme(plot.title = element_text(hjust = 0.5))+
							coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE)
						legend.p <- NULL
						if(vector.friendly){
							####p <- Seurat::AugmentPlot(p,width=width,height=height)
							tmpfilename <- sprintf("%s.tmp.%s.png",if(!is.null(out.prefix)) out.prefix else "",x)
							ggsave(filename=tmpfilename, plot = p + theme_void() + theme(legend.position = "none",
												  axis.line.x = element_blank(), axis.line.y = element_blank(),
												  axis.title.x = element_blank(), axis.title.y = element_blank(),
												  axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
												  plot.title = element_blank()),
								 width=7,height=6)
							pbuild.params <- ggplot_build(plot = p)$layout$panel_params[[1]]
							range.values <- c( pbuild.params$x.range, pbuild.params$y.range)
							img <- png::readPNG(source = tmpfilename)
							##print(str(p$data))
							##p.test <<- p
							blank <- ggplot(data = as.data.table(p$data)[1,,drop=F],aes(Dim1,Dim2)) +
									  geom_blank()+
									  labs(title=x, x ="", y = "") +
									  theme(plot.title = element_text(hjust = 0.5))
							blank <- blank + p$theme + coord_cartesian(xlim = range.values[1:2], ylim = range.values[3:4], expand = F)
							blank <- blank + annotation_raster(raster = img,
															 xmin = range.values[1], xmax = range.values[2],
															 ymin = range.values[3], ymax = range.values[4])
						    legend.blank <- cowplot::get_legend(p)
							if(show.legend){
								p <- cowplot::plot_grid(blank, legend.blank, rel_widths = c(3, 1))
							}else{
								p <- blank
							}
							legend.p <- legend.blank
							file.remove(tmpfilename)
						}else{
							legend.p <- cowplot::get_legend(p)
							if(!show.legend){
								p <- p + theme(legend.position = "none")
							}
						}
						return(list("plot"=p,"legend"=legend.p))
					 })
  }

  if(scales=="fixed"){
      if(!is.null(clamp) && clamp=="none"){
      }else{
          if(is.null(clamp)){
              clamp <- quantile(dat.plot.melt$value,c(0.05,0.95))
          }
          dat.plot.melt[value < clamp[1],value:=clamp[1]]
          dat.plot.melt[value > clamp[2],value:=clamp[2]]
      }
	  value.range <- pretty(dat.plot.melt$value)
      multi.p <- make.plot(gene.to.show,dat.plot.melt,show.legend=F,clamp="none",value.range=value.range,size=size,vector.friendly=vector.friendly,out.prefix=out.prefix)
	  legend.p <- multi.p[[1]][["legend"]]
	  ####+ theme(legend.box.margin = margin(0, 0, 0, 12)))
      p <- cowplot::plot_grid(plotlist=llply(multi.p,function(x){ x[["plot"]] }),ncol = p.ncol,align = "hv")
	  p <- cowplot::plot_grid(p, legend.p, rel_widths = c(3, .4))
  }else{
      multi.p <- make.plot(gene.to.show,dat.plot.melt,show.legend=T,clamp=clamp,value.range=NULL,size=size,vector.friendly=vector.friendly,out.prefix=out.prefix)
      p <- cowplot::plot_grid(plotlist=llply(multi.p,function(x){ x[["plot"]] }),ncol = p.ncol,align = "hv")
  }
  if(!is.null(out.prefix)){
    ggplot2::ggsave(sprintf("%s.geneOntSNE.pdf",out.prefix),plot=p,width = width,height = height)
  }
  return(p)
}

#' Wrap for plotting 2D density
#'
#' @importFrom ks kde
#' @importFrom fields image.plot
#' @param x matrix or data.frame; map data, row for sample, column for dimension
#' @param peaks integer or character; index or names of the peaks. (default: NULL)
#' @usage plot.density2D(x, peaks)
#' @details use ks::kde for density estimation
#'
plot.density2D <- function(x,peaks=NULL)
{
  .density <- ks::kde(x)
  ##dev.new()
  par(mar=c(5,4,5,6))
  .zz <- c(10,20,30,40,50,60,70,80,90)
  plot(.density,display="filled.contour2", cont=.zz,xlab="Dim1", ylab="Dim2")
  fields::image.plot(zlim=c(0,.zz[length(.zz)]),legend.only=TRUE, col = c("transparent", rev(heat.colors(length(.zz)))),
                     axis.args=list( at=.zz, labels=sprintf("%s%%",100-.zz)), legend.width=2.0,legend.mar=4.5)
  if(!is.null(peaks)){
    plot(.density,display="filled.contour2", cont=.zz,xlab="Dim1", ylab="Dim2")
    points(x[peaks,,drop=F],pch=3,cex=2,col="black")
    fields::image.plot(zlim=c(0,.zz[length(.zz)]),legend.only=TRUE, col = c("transparent", rev(heat.colors(length(.zz)))),
                       axis.args=list( at=.zz, labels=sprintf("%s%%",100-.zz)), legend.width=2.0,legend.mar=4.5)
  }
  ##pp <- recordPlot()
  ##dev.off()
  #pp
}


#' plot matrix (typically genes expression)
#' @param dat matrix; matrix
#' @param out.prefix character; output prefix.
#' @param mytitle character; (default: "Heatmap")
#' @param show.number logical; (default: NULL)
#' @param do.clust logical, character or dendrogram; passed to both cluster_columns and cluster_rows of Heatmap. Higher priority than clust.row and clust.column (default: NULL)
#' @param clust.row logical, character or dendrogram; passed to cluster_rows of Heatmap (default: FALSE)
#' @param clust.column logical, character or dendrogram; passed to cluster_columns of Heatmap (default: FALSE)
#' @param waterfall.row logical, order rows to make plot like waterfall (default: FALSE)
#' @param waterfall.column logical, order rows to make plot like waterfall (default: FALSE)
#' @param show.dendrogram logical, whetehr show the dendrogram (default: FALSE)
#' @param z.lo double; (default: NULL)
#' @param z.hi double; (default: NULL)
#' @param z.len integer; (default: 100)
#' @param col.ht vector; (default: NULL)
#' @param palatte character; (default: NULL)
#' @param row.ann.dat data.frame; data for row annotation; (default: NULL)
#' @param row.split vector; used for row; (default: NULL)
#' @param returnHT logical; whether return HT; (default: FALSE)
#' @param par.legend list; lengend parameters, used to overwrite the default setting; (default: list())
#' @param par.heatmap list; other heatmap parameters, (default: list())
#' @param par.warterfall list; parameters for warterfall, sucah as score.alpha
#' @param pdf.width double; width of the output plot (default: 22)
#' @param pdf.height double; height of the output plot (default: 22)
#' @param fig.type character; type of output file (default: "pdf")
#' @param exp.name character; showd in the legend (default: "Count")
#' @param ... parameters passed to ComplexHeatmap::draw
#' @import data.table
#' @import ggplot2
#' @importFrom ggpubr ggscatter
#' @importFrom stats dist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @details plot matrix
#' @export
plot.matrix.simple <- function(dat,out.prefix=NULL,mytitle="Heatmap",show.number=NULL,
                               do.clust=NULL,z.lo=NULL,z.hi=NULL,z.len=100,palatte=NULL,
                               clust.row=FALSE,clust.column=FALSE,show.dendrogram=FALSE,
                               waterfall.row=FALSE,waterfall.column=FALSE,
                               row.ann.dat=NULL,row.split=NULL,returnHT=FALSE,
                               par.legend=list(),par.heatmap=list(),col.ht=NULL,
							   par.warterfall=list(score.alpha=1.5,do.norm=T),
                               pdf.width=8,pdf.height=8,fig.type="pdf",exp.name="Count",...)
{
    require("gplots")
    require("ComplexHeatmap")
    require("circlize")
    require("gridBase")
    require("RColorBrewer")

    #### ordering the matrix
    if(!is.null(do.clust)){
        clust.row <- do.clust
        clust.column <- do.clust
    }
    if(!is.null(row.split) && is.null(names(row.split))){
        names(row.split) <- rownames(dat)
    }
	dat.list <- do.call(matrix.waterfall,c(list(dat=dat,clust.row=clust.row,clust.column=clust.column,
												waterfall.column=waterfall.column,waterfall.row=waterfall.row,
												type.return="list"),par.warterfall))
	dat <- dat.list[["dat"]]
	clust.column <- dat.list[["clust.column"]]
	clust.row <- dat.list[["clust.row"]]
	scoresR <- dat.list[["scoresR"]]
	scoresC <- dat.list[["scoresC"]]
    if(!is.null(row.ann.dat)){ row.ann.dat <- row.ann.dat[order(scoresR,decreasing = T),] }

	if(!is.null(out.prefix)){
		if(fig.type=="pdf"){
			pdf(sprintf("%s.pdf",out.prefix),width=pdf.width,height=pdf.height)
		}else{
			png(sprintf("%s.png",out.prefix),width=pdf.width*100,height=pdf.height*100)
		}
        opar <- par(mar=c(4,2,4,4))
        plot.new()
        title(main = mytitle,cex.main=2)
        #legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)
        ### Integrating Grid Graphics Output with Base Graphics Output
        vps <- baseViewports()
        pushViewport(vps$inner, vps$figure, vps$plot)
    }
    tmp.var <- pretty((dat),n=8)
    if(is.null(z.lo)){ z.lo <- tmp.var[1] }
    if(is.null(z.hi)){ z.hi <- tmp.var[length(tmp.var)] }

    my.cell_fun <- NULL
	if(!is.null(show.number)){
		if(is.logical(show.number) && show.number==T){
			my.cell_fun <- function(j, i, x, y, w, h, col) { grid.text(dat[i, j], x, y) }
		}else if(is.matrix(show.number)){
			my.cell_fun <- function(j, i, x, y, w, h, col) { grid.text(show.number[i, j], x, y) }
		}
	}
    m <- ncol(dat)
    n <- nrow(dat)
    if(is.null(palatte)){
        palatte <- rev(brewer.pal(n = 7,name = "RdYlBu"))
    }
    par.legend.used <- list(title = exp.name,
                               grid_width = unit(0.4, "cm"),
                               grid_height = unit(0.4, "cm"),
                               #legend_width=2,
                               legend_height=unit(6,"cm"),
                               title_gp = gpar(fontsize = 10, fontface = "bold"),
                               #color_bar = "continuous",
                               label_gp = gpar(fontsize = 8))
    if(length(par.legend)>0){
        for(ipar in names(par.legend)){
            par.legend.used[[ipar]] <- par.legend[[ipar]]
        }
    }
	par.heatmap.used <- list(cex.row=1,cex.column=1)
	if(length(par.heatmap)>0){
		for(ipar in names(par.heatmap)){
			par.heatmap.used[[ipar]] <- par.heatmap[[ipar]]
		}
	}
	.cex.row <- par.heatmap.used[["cex.row"]]
	.cex.column <- par.heatmap.used[["cex.column"]]

    if(!is.null(row.split)){
        row.split <- row.split[rownames(dat)]
    }

    ht <- ComplexHeatmap::Heatmap(dat, name = mytitle,
                  col = if(is.null(col.ht)) colorRamp2(seq(z.lo,z.hi,length=z.len), colorRampPalette(palatte)(z.len)) else col.ht,
                  cluster_columns=clust.column,cluster_rows=clust.row,
                  row_dend_reorder = FALSE, column_dend_reorder = FALSE,
                  column_names_gp = gpar(fontsize = 12*28*.cex.column/max(m,32)),
                  row_names_gp = gpar(fontsize = 10*28*.cex.row/max(n,32)),
                  #row_dend_width = unit(4, "cm"),
                  #column_dend_height = unit(4, "cm"),
                  show_row_dend = show.dendrogram,
                  show_column_dend = show.dendrogram,
                  heatmap_legend_param = par.legend.used,
                  cell_fun = my.cell_fun,...)
    if(!is.null(row.ann.dat)){
        for(idx in colnames(row.ann.dat)){
            idx.col <- NULL
            if(is.logical(row.ann.dat[[idx]])){
                row.ann.dat[[idx]] <- as.character(row.ann.dat[[idx]])
                idx.col <- c("TRUE" = "red", "FALSE" = "blue")
            }else if(is.character(row.ann.dat[[idx]])){
                idx.levels <- sort(unique(row.ann.dat[[idx]]))
                idx.col <- structure(sscClust:::auto.colSet(length(idx.levels)),names=idx.levels)
            }
            vv <- row.ann.dat[[idx]]
            names(vv) <- rownames(dat)
            ht <- ht + ComplexHeatmap::Heatmap(vv,name=idx,col=idx.col,
                                               row_names_gp = gpar(fontsize = 10*28/max(n,32)))
        }
    }

    if(!is.null(out.prefix)){
        ComplexHeatmap::draw(ht, newpage= FALSE,merge_legends = TRUE,split=row.split)
        dev.off()
        #par(opar)
    }
    if(returnHT){ return(ht) }
}

#' plot branch
#' @param obj.clust object;
#' @param out.prefix character; output prefix.
#' @param ncls integer; (default: 1)
#' @param cluster integer vector; (default: NULL)
#' @importFrom  dendextend color_branches set
#' @importFrom moduleColor plotHclustColors
#' @details plot dendrogram
plot.branch <- function(obj.clust,out.prefix,ncls=1,cluster=NULL)
{
    obj.dend <- as.dendrogram(obj.clust)
    if(!is.null(cluster)){
        ncls <- length(unique(cluster))
        colSet.cls <- auto.colSet(ncls,"Paired")
        names(colSet.cls) <- unique(cluster[order.dendrogram(obj.dend)])
        col.cls <- data.frame("k0"=sapply(cluster,function(x){ colSet.cls[as.character(x)] }),
                                stringsAsFactors=F)
        branch.col <- color_branches(obj.dend,
                                     clusters=cluster[order.dendrogram(obj.dend)],
                                     col=colSet.cls)
    }else{
        dend.cutree <- cutree(obj.clust, c(ncls,ncls), order_clusters_as_data = T)
        colSet.cls <- auto.colSet(ncls,"Paired")
        col.cls <- t(apply(dend.cutree,1,function(x){ colSet.cls[x] }))
        branch.col <- color_branches(obj.dend,k=ncls,col=colSet.cls)
        colnames(col.cls) <- c("k0","k0")
        col.cls <- col.cls[,1,drop=F]
    }
    
    branch.col <- branch.col %>%
        dendextend::set("branches_lwd", 1.5) %>%
        dendextend::set("labels_colors",col.cls$k0[order.dendrogram(obj.dend)]) %>%
        dendextend::set("labels_cex", 1*50/max(length(obj.clust$labels),32))

    pdf(sprintf("%s.branch.pdf",out.prefix),width=12,height=8)
    layout(matrix(c(1,2),nrow = 2),heights = c(0.8,0.2))
    par(mar=c(15,4,4,2),xpd=T)
    plot(branch.col)
    #par(mar=c(0,4,4,2),xpd=T)
    #plot(obj.clust,sub="",xlab="",hang=-1,cex=1.0*50/max(length(obj.clust$labels),32))
    par(mar=c(5,4,0,2))
    moduleColor::plotHclustColors(obj.clust, colors=col.cls, cex.rowLabels = 1.1)
    dev.off()

    return(list("obj.clust"=obj.clust,"branch"=branch.col))
}

#' plot distribution using cell info table
#' @param obj object; can be class of Seurat, SingleCellExperiment or data.frame;
#' @param out.prefix character; output prefix.
#' @param plot.type character; (default: "barplot")
#' @param facet.ncol integer; (default: 3)
#' @param plot.width integer; (default: 10)
#' @param plot.height integer; (default: 5)
#' @param test.method character; (default: "fisher.test")
#' @param group.filter character; (default: NULL)
#' @param sort.freq logical; (default: FALSE)
#' @param bar.position ; passed to ggbarplot's position (default: position_dodge2())
#' @param cmp.var character; (default: "Species")
#' @param group.var character; (default: "ClusterID")
#' @param donor.var character; (default: "donor")
#' @param verbose logical; (default: FALSE)
#' @param par.stat parameter passed to stat_compare_means (default: list() )
#' @param min.NTotal minimum number of cells to calculate frequency (default: 0 )
#' @param ... parameter passed to ggpubr
#' @importFrom  ggpubr ggbarplot ggboxplot stat_compare_means
#' @importFrom ggplot2 facet_wrap coord_cartesian expand_limits geom_text element_text theme ggsave
#' @importFrom plyr ldply
#' @details plot distribution
#' @export
plotDistFromCellInfoTable <- function(obj,out.prefix,plot.type="barplot",
                        facet.ncol=3,plot.width=10,plot.height=5,test.method="fisher.test",
						group.filter=NULL,sort.freq=F,bar.position=position_dodge2(),
						par.stat=list(),min.NTotal=0,
                        cmp.var="Species",group.var="ClusterID",donor.var="donor",verbose=F,...)
{
    require("ggpubr")
    require("ggplot2")
    require("plyr")
    require("ggsignif")
    if("Seurat" %in% class(obj)){
        dat.tb <- as.data.table(obj[[]])
    }else if("SingleCellExperiment" %in% class(obj)){
        dat.tb <- as.data.table(colData(obj))
    }else if(is.data.frame(obj)){
        dat.tb <- as.data.table(obj)
    }
    if(!is.factor(dat.tb[[cmp.var]])){
        dat.tb[[cmp.var]] <- factor(dat.tb[[cmp.var]])
    }
    if(is.null(donor.var)){
        dat.spe.group.dist <- dat.tb[, .N, by=c(cmp.var,group.var)]
        colnames(dat.spe.group.dist) <- c("cmp.var","group.var","N")
        dat.spe.group.dist$donor.var <- "ALL"
    }else{
        dat.spe.group.dist <- dat.tb[, .N, by=c(donor.var,cmp.var,group.var)]
        colnames(dat.spe.group.dist) <- c("donor.var","cmp.var","group.var","N")
    }
    #### fill missing value using 0
    dat.spe.group.dist[,donor.var:=gsub("_",".",donor.var)]
    dat.spe.group.dist[,cmp.var:=gsub("_",".",cmp.var)]
    dat.spe.group.dist <- melt(dcast(dat.spe.group.dist,group.var~donor.var+cmp.var,fill=0,value.var="N"),
                               id.vars="group.var",value.name="N")
    dat.spe.group.dist[,donor.var:=sapply(strsplit(as.character(variable),"_",perl=T),"[",1) ]
    dat.spe.group.dist[,cmp.var:=sapply(strsplit(as.character(variable),"_",perl=T),"[",2) ]
    dat.spe.group.dist[,variable:=NULL]
    dat.spe.group.dist <- dat.spe.group.dist[,.(group.var=group.var,N=N,NTotal=sum(.SD$N)),
                                             by=c("donor.var","cmp.var")]
    dat.spe.group.dist[,freq:=N/NTotal]
	##### data for plot: dat.spe.group.dist
	dat.spe.group.dist <- dat.spe.group.dist[NTotal>=min.NTotal,]
	###dat.spe.group.dist.test <<- dat.spe.group.dist
	dat.freq.med <- dat.spe.group.dist[,.(freq.med=median(freq)),by=c("group.var","cmp.var")][order(freq.med),]
	if(!is.null(group.filter)){
		dat.spe.group.dist <- dat.spe.group.dist[group.var==group.filter,]
		dat.freq.med <- dat.freq.med[group.var==group.filter]
	}
	if(sort.freq){
		dat.spe.group.dist <- dat.spe.group.dist[order(freq,cmp.var),]
		dat.spe.group.dist[,cmp.var:=factor(cmp.var,levels=unique(dat.freq.med$cmp.var))]
	}else if(is.factor(dat.tb[[cmp.var]])){
        clevels <- levels(dat.tb[[cmp.var]])
        dat.spe.group.dist[,cmp.var:=factor(cmp.var,levels=clevels)]
	}

    if(plot.type=="boxplot"){
        p <- ggboxplot(dat.spe.group.dist,x="group.var",y="freq",
                       color = "cmp.var", 
                       add = "jitter",outlier.shape=NA,...)
        if(test.method!=""){
			p <- p + do.call(stat_compare_means,c(list(mapping=aes_string(group = "cmp.var"),
													   label="p.signif"),par.stat))
        }
    }else if(plot.type=="boxplot2"){
        p <- ggboxplot(dat.spe.group.dist,x="cmp.var",y="freq",
                       color = "cmp.var", 
                       add = "jitter",outlier.shape=NA,...) +
                facet_wrap(~group.var,ncol=facet.ncol,scales="free_y")+
                expand_limits(y=0)
        if(test.method!=""){
			p <- p + do.call(stat_compare_means,c(list(label="p.format"),par.stat))
        }
    }else if(plot.type=="barplot"){
        dat.plot <- dat.spe.group.dist[,.(N=sum(.SD$N),NTotal=sum(.SD$NTotal),
                                          freq=mean(.SD$freq)),
                                        by=c("cmp.var","group.var")]
        dat.plot[,NOther:=NTotal-N]
        
        p <- ggbarplot(dat.plot,x="group.var",y="freq",fill="cmp.var",color=NA,
                       position=bar.position,...)

        if(test.method!=""){
            if(test.method=="prop.test"){
                ann.tb <- dat.plot[,.(p.value=prop.test(.SD$N, .SD$NTotal)$p.value,
                                      y_pos=max(.SD$freq)),
                                   by="group.var"]
            }else if(test.method=="fisher.test"){
                ann.tb <- dat.plot[,.(p.value=fisher.test(.SD[,c("N","NOther"),with=F])$p.value,
                                      y_pos=max(.SD$freq)),
                                   by="group.var"]
            }
            ann.tb[,p.adj:=p.adjust(p.value,"fdr")]
            ann.tb[,p.signif:="ns"]
            ann.tb[p.adj<0.05,p.signif:="*"]
            ann.tb[p.adj<0.01,p.signif:="**"]
            ann.tb[p.adj<0.001,p.signif:="***"]
            #ann.tb[,xmin:=as.numeric(factor(group.var))-0.2]
            #ann.tb[,xmax:=as.numeric(factor(group.var))+0.2]

            p <- p + geom_text(data=ann.tb,aes(x=group.var,y=y_pos,label=p.signif),vjust=-0.5)
        }

    }

	if(plot.type=="none"){
		return(dat.spe.group.dist)
	}else{
		p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
					coord_cartesian(clip="off")
		if(verbose){
			write.table(dat.spe.group.dist, sprintf("%s.dist.%s.%s.freq.txt",out.prefix,cmp.var,group.var),
						row.names=F,sep="\t",quote=F)
		}
		if(!is.null(out.prefix)){
			ggsave(sprintf("%s.dist.%s.%s.%s.freq.pdf",out.prefix,plot.type,cmp.var,group.var),
				   width=plot.width,height=plot.height)
		}else{
			return(p)
		}
	}
}

#' make the matrix looks like "waterfall" (typically genes expression)
#' @param dat matrix; matrix
#' @param clust.row logical, character or dendrogram; passed to cluster_rows of Heatmap (default: FALSE)
#' @param clust.column logical, character or dendrogram; passed to cluster_columns of Heatmap (default: FALSE)
#' @param waterfall.row logical, order rows to make plot like waterfall (default: FALSE)
#' @param waterfall.column logical, order rows to make plot like waterfall (default: FALSE)
#' @param score.alpha double; for row/column score; (default: 1.5)
#' @param do.norm logical; normalized the data before calculating scores?; (default: TRUE)
#' @param type.return character; return type; (default: "matrix")
#' @param ... parameters passed to run.cutreeDynamic
#' @import data.table
#' @import ggplot2
#' @importFrom stats dist
#' @export
matrix.waterfall <- function(dat,score.alpha=1.5,clust.row=FALSE,clust.column=FALSE,
							 waterfall.column=F,waterfall.row=F,
							 type.return="matrix",do.norm=T,...)
{
    scoreVec = function(x) {
        score = 0
        if(all(x<1)){ x <- x^100 }
		if(do.norm){ x <- x/sum(x) }
        m <- length(x)
        score <- sum(score.alpha^(-seq_len(m)) * x)
        if(is.na(score)) { score <- 0 }
        return(score)
    }
	scoresC <- NULL
	scoresR <- NULL
    dat.ordered <- dat
    if(clust.row=="cutreeDynamic"){
        res.clust.row <- run.cutreeDynamic(dat,...)
        dat.ordered <- dat[res.clust.row$hclust$order,]
        clust.row <- res.clust.row$branch
    }else if(clust.row=="cutree"){
        res.clust.row <- run.cutree(dat,...)
        dat.ordered <- dat[res.clust.row$hclust$order,]
        clust.row <- res.clust.row$branch
	}
    if(clust.column=="cutreeDynamic"){
        res.clust.column <- run.cutreeDynamic(t(dat),...)
        dat.ordered <- dat[,res.clust.column$hclust$order]
        clust.column <- res.clust.column$branch
    }else if(clust.column=="cutree"){
        res.clust.column <- run.cutree(t(dat),...)
        dat.ordered <- dat[,res.clust.column$hclust$order]
        clust.column <- res.clust.column$branch
	}
    if(waterfall.column){
        scoresC <- apply(dat.ordered, 2, scoreVec)
        dat <- dat[,order(scoresC,decreasing = T)]
    }
    if(waterfall.row){
        scoresR <- apply(dat.ordered, 1, scoreVec)
        dat <- dat[order(scoresR,decreasing = T),]
    }
	if(type.return=="matrix"){
		return(dat)
	}else if(type.return=="list"){
		return(list("dat"=dat,
					"clust.column"=clust.column,"clust.row"=clust.row,
					"scoresC"=scoresC,"scoresR"=scoresR))
	}
}



