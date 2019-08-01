#!/usr/bin/env Rscript

library("Seurat")
library("magrittr")
library("tictoc")
library("plyr")
library("tibble")
library("dplyr")
library("doParallel")
library("sscClust")
library("limma")
library("data.table")
library("cowplot")
library("ggplot2")
library("ggpubr")
library("R.utils")
library("data.table")
library("cowplot")
library("dendextend")
library("moduleColor")

exp.human.file <- "baron.sce.flt.hom.human.rds"
exp.mouse.file <- "baron.sce.flt.hom.mouse.rds"
deg.human.file <- "OUT.cor/baron.cor.deg.human.txt"
deg.mouse.file <- "OUT.cor/baron.cor.deg.mouse.txt"

out.prefix <- "./OUT.seurat3/baron.seurat"
dir.create(dirname(out.prefix),F,T)

################

#if(F){
sce.human <- readRDS(exp.human.file)
sce.mouse <- readRDS(exp.mouse.file)

#gene.de.list <- list()
#gene.de.list[["human"]] <- fread(deg.human.file)
#gene.de.list[["mouse"]] <- fread(deg.mouse.file)
#gene.de.list[["human"]][1:4,1:6]
#gene.de.list[["mouse"]][1:4,1:6]
#gene.de.union.vec <- unique(c(gene.de.list[["human"]]$geneID,
#                              gene.de.list[["mouse"]]$geneID))
#gene.common <- intersect(gid.mapping[rownames(sce.human)],gene.de.union.vec)
#gene.common <- intersect(gene.common,rownames(sce.mouse))

gene.common <- intersect(rownames(sce.human),rownames(sce.mouse))

seu.list <- list()

tmp.dat <- assay(sce.human,"logcounts")[gene.common,]
tmp.dat[1:4,1:5]
seu.list[["human"]] <- CreateSeuratObject(counts=tmp.dat,
                                          meta.data=as.data.frame(colData(sce.human)))
seu.list[["human"]] <- SetAssayData(seu.list[["human"]], slot = "data", new.data = tmp.dat)
seu.list$human$majorCluster <- sprintf("human.%s",seu.list$human$cell_type1)
table(seu.list$human$majorCluster)
seu.list[["human"]]

tmp.dat <- assay(sce.mouse,"logcounts")[gene.common,]
tmp.dat[1:4,1:5]
seu.list[["mouse"]] <- CreateSeuratObject(counts=tmp.dat,
                                          meta.data=as.data.frame(colData(sce.mouse)))
seu.list[["mouse"]] <- SetAssayData(seu.list[["mouse"]], slot = "data", new.data = tmp.dat)
seu.list$mouse$majorCluster <- sprintf("mouse.%s",seu.list$mouse$cell_type1)
table(seu.list$mouse$majorCluster)
seu.list[["mouse"]]

all(rownames(seu.list$human)==rownames(seu.list$mouse))

for (i in 1:length(x = seu.list)) {
###    seu.list[[i]] <- NormalizeData(object = seu.list[[i]], verbose = FALSE)
    seu.list[[i]] <- FindVariableFeatures(object = seu.list[[i]], selection.method = "vst",
                                          nfeatures = 2000, verbose = FALSE)
}


####### M1. data integration

tic("FindIntegrationAnchors")
seu.anchors <- FindIntegrationAnchors(object.list = seu.list, dims = 1:30)
toc()

tic("IntegrateData")
seu.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:30)
toc()

DefaultAssay(object = seu.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
seu.integrated <- ScaleData(object = seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(object = seu.integrated, npcs = 30, verbose = FALSE)
seu.integrated <- RunUMAP(object = seu.integrated, reduction = "pca", dims = 1:30)

m <- regexec("^(.+)\\.",seu.integrated$majorCluster,perl=T)
mm <- regmatches(seu.integrated$majorCluster,m)
seu.integrated$species <- sapply(mm,"[",2)

saveRDS(seu.integrated,sprintf("%s.integrated.seu.rds",out.prefix))
saveRDS(seu.anchors,sprintf("%s.integrated.anchors.rds",out.prefix))

#### side by side comparison


p1 <- DimPlot(object = seu.integrated,
              reduction = "umap",
              group.by = "species",label=TRUE)
ggsave(sprintf("%s.umap.species.pdf",out.prefix),width=5,height=4)
p1 <- DimPlot(object = seu.integrated,
              reduction = "umap",
              group.by = "species", split.by="species",label=TRUE)
ggsave(sprintf("%s.umap.species.split.pdf",out.prefix),width=8,height=4)

p1 <- DimPlot(object = seu.integrated,
              reduction = "umap",
              group.by = "majorCluster",label=F)
ggsave(sprintf("%s.umap.majorCluster.pdf",out.prefix),width=9,height=4)

p1 <- DimPlot(object = seu.integrated[,seu.integrated$species=="human"],
              reduction = "umap",
              group.by = "majorCluster",label=T,label.size=2,repel=T)
ggsave(sprintf("%s.umap.human.majorCluster.pdf",out.prefix),width=7,height=4)

p1 <- DimPlot(object = seu.integrated[,seu.integrated$species=="mouse"],
              reduction = "umap",
              group.by = "majorCluster",label=T,label.size=2,repel=T)
ggsave(sprintf("%s.umap.mouse.majorCluster.pdf",out.prefix),width=7,height=4)

#######################################################

#### clustring
seu.integrated <- FindNeighbors(object = seu.integrated, reduction = "pca", dims = 1:30)

resolution.vec <- c(0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2)
seu.integrated <- FindClusters(object = seu.integrated,resolution = resolution.vec)

for(opt.res in resolution.vec){
    print(table(seu.integrated[[sprintf("integrated_snn_res.%s",opt.res)]]))
}

### plot
plot.resolution.list <- list()
for(opt.res in resolution.vec){
    plot.resolution.list[[sprintf("integrated_snn_res.%s",opt.res)]] <- DimPlot(seu.integrated, reduction = "umap",
                                                                         pt.size=0.15,label.size=2,
                                                                         label = TRUE,
                                                                         group.by=sprintf("integrated_snn_res.%s",
                                                                                              opt.res)) +
        NoLegend()
}

pp <- plot_grid(plotlist=plot.resolution.list[1:4],
                ncol = 2,align = "hv")
save_plot(sprintf("%s.umap.res.1.png",out.prefix),pp, ncol = 2, base_aspect_ratio=0.55)
pp <- plot_grid(plotlist=plot.resolution.list[5:8],
                ncol = 2,align = "hv")
save_plot(sprintf("%s.umap.res.2.png",out.prefix),pp, ncol = 2, base_aspect_ratio=0.55)

### optimal res: 
### 1.0
seu.integrated$ClusterID <- sprintf("C%02d",as.integer(as.character(seu.integrated$integrated_snn_res.1)))

write.table(cbind(data.frame(sample=rownames(seu.integrated@meta.data),stringsAsFactors = F),
                  seu.integrated@meta.data),
            sprintf("%s.clusterInfo.txt",out.prefix),sep = "\t",row.names = F,quote = F)
saveRDS(seu.integrated,sprintf("%s.integrated.rds",out.prefix))
###seu.integrated <- readRDS(sprintf("%s.integrated.rds",out.prefix))

p <- DimPlot(object = seu.integrated, reduction = "umap", group.by = "ClusterID",label=TRUE)
ggsave(sprintf("%s.umap.ClusterID.pdf",out.prefix),width=6,height=5)



