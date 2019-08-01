#!/usr/bin/env Rscript

library("sscClust")
library("plyr")
library("tictoc")
library("data.table")

exp.human.file <- "baron.sce.flt.hom.human.rds"
exp.mouse.file <- "baron.sce.flt.hom.mouse.rds"
out.prefix <- "./OUT.cor/baron.cor"
dir.create(dirname(out.prefix),F,T)

sce.list <- list()
sce.list[["human"]] <- readRDS(exp.human.file)
sce.list[["mouse"]] <- readRDS(exp.mouse.file)
#sce.list$human$dataset.id <- "human"
#sce.list$mouse$dataset.id <- "mouse"

tic("integrate.by.avg")
sce.pb <- integrate.by.avg(sce.list,out.prefix,assay.name="logcounts",ncores=6,avg.by="cell_type1")
toc()

comb.gene.table <- rbind(as.data.table(metadata(sce.pb)$ssc$gene.de.list$human[,1:8]),
                         as.data.table(metadata(sce.pb)$ssc$gene.de.list$mouse[,1:8]))
comb.gene.table <- comb.gene.table[order(cluster,-AUC),]
ff <- !duplicated(comb.gene.table$geneID)
comb.gene.table.flt <- comb.gene.table[ff,][,head(.SD,n=5),by=c("cluster")]
comb.gene.table.flt$Group <- comb.gene.table.flt$cluster

obj.hclust <-  metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$hclust
comb.gene.table.flt$Group <- factor(comb.gene.table.flt$Group,
                                    levels=unique(gsub("^.+\\.","",obj.hclust$labels[obj.hclust$order])))
comb.gene.table.flt <- comb.gene.table.flt[order(Group),]

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 6)

ssc.plot.heatmap(sce.pb,out.prefix=sprintf("%s.exp.avg.example.preOrderCol",out.prefix),
                 gene.desc=comb.gene.table.flt,
                 columns="pca.dynamicTreeCut.kauto",
                 do.scale=F,z.hi=1.5,z.lo=-1.5,
                 do.clustering.col=T,
                 do.clustering.row=F,
                 pdf.width=14,pdf.height=16,
                 dend.col=metadata(sce.pb)$ssc$clust.res$dynamicTreeCut$auto$branch
                 )

saveRDS(sce.pb,sprintf("%s.sce.pb.rds",out.prefix))
#sce.pb <- readRDS(sprintf("%s.sce.pb.rds",out.prefix))

write.table(metadata(sce.pb)$ssc$gene.de.list$human,sprintf("%s.deg.human.txt",out.prefix),
            sep="\t",quote=F,row.names=F)
write.table(metadata(sce.pb)$ssc$gene.de.list$mouse,sprintf("%s.deg.mouse.txt",out.prefix),
            sep="\t",quote=F,row.names=F)


