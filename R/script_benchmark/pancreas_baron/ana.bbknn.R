#!/usr/bin/env Rscript

library("reticulate")
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

out.prefix <- "./OUT.bbknn/baron.bbknn"
dir.create(dirname(out.prefix),F,T)

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
#gene.common <- intersect(gid.mapping[rownames(sce.human)],rownames(sce.mouse))

gene.common <- intersect(rownames(sce.human),rownames(sce.mouse))

seu.list <- list()

tmp.dat <- assay(sce.human,"logcounts")[gene.common,]
tmp.dat[1:4,1:5]
seu.list[["human"]] <- CreateSeuratObject(counts=tmp.dat,
                                          meta.data=as.data.frame(colData(sce.human)))
seu.list[["human"]] <- SetAssayData(seu.list[["human"]], slot = "data", new.data = tmp.dat)
seu.list$human$majorCluster <- sprintf("human.%s",seu.list$human$cell_type1)
table(seu.list$human$majorCluster)
seu.list[["human"]]$batch <- "human"
seu.list[["human"]]


tmp.dat <- assay(sce.mouse,"logcounts")[gene.common,]
tmp.dat[1:4,1:5]
seu.list[["mouse"]] <- CreateSeuratObject(counts=tmp.dat,
                                          meta.data=as.data.frame(colData(sce.mouse)))
seu.list[["mouse"]] <- SetAssayData(seu.list[["mouse"]], slot = "data", new.data = tmp.dat)
seu.list$mouse$majorCluster <- sprintf("mouse.%s",seu.list$mouse$cell_type1)
table(seu.list$mouse$majorCluster)
seu.list[["mouse"]]$batch <- "mouse"
seu.list[["mouse"]]

all(rownames(seu.list$human)==rownames(seu.list$mouse))

for (i in 1:length(x = seu.list)) {
###    seu.list[[i]] <- NormalizeData(object = seu.list[[i]], verbose = FALSE)
    seu.list[[i]] <- FindVariableFeatures(object = seu.list[[i]], selection.method = "vst",
                                          nfeatures = 2000, verbose = FALSE)
}

gene.use <- unique(c(VariableFeatures(seu.list$human),
                     VariableFeatures(seu.list$mouse)))

#### reticulate
anndata = import("anndata",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
sc = import("scanpy",convert=FALSE)
plt = import("matplotlib.pyplot",convert=FALSE)
plt$switch_backend('Agg')

dat.py.human <- anndata$AnnData(X=t(GetAssayData(seu.list$human,"data")[gene.use,]),
                                obs=seu.list$human[[]])
dat.py.human$var_names <- gene.use

dat.py.mouse <- anndata$AnnData(X=t(GetAssayData(seu.list$mouse,"data")[gene.use,]),
                                obs=seu.list$mouse[[]])
dat.py.mouse$var_names <- gene.use

adata = dat.py.human$concatenate(dat.py.mouse, join='outer')

##sc.pp.log1p(adata)
sc$pp$scale(adata, max_value=10)
sc$tl$pca(adata)
#adata$obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc$pl$pca_variance_ratio(adata, log=TRUE)

num_pcs = 20L

adata_bbknn_trim = bbknn$bbknn(adata, neighbors_within_batch=5L, n_pcs=num_pcs, trim=50L, copy=TRUE)
sc$tl$umap(adata_bbknn_trim)
sc$tl$louvain(adata_bbknn_trim)

adata_bbknn_trim$write(sprintf("%s.bbknn_trim.h5ad",out.prefix))
adata_bbknn_trim$obs$to_csv(sprintf("%s.bbknn_trim.obs.csv",out.prefix))

#umap = py_to_r(adata_bbknn_trim$obsm$X_umap)

sc$pl$umap(adata_bbknn_trim, color=c('batch','louvain','majorCluster'))
fig = plt$gcf()
fig$set_size_inches(15, 4)
#fig$tight_layout()
fig$savefig(sprintf("%s.umap.species.pdf",out.prefix), dpi=100,bbox_inches="tight",pad_inches=0.5)




