#!/usr/bin/env Rscript

library("sscClust")
library("plyr")


geneID.homolog <- function(gid,gid.type="entrezgene")
{
    require("data.table")
    require("biomaRt")
    
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

    g.list.human.ext = as.data.table(getBM(attributes = c('ensembl_gene_id',
                                            "entrezgene",
                                            "hgnc_symbol",
                                            "gene_biotype",
                                            "description",
                                            "chromosome_name",
                                            "start_position",
                                            "end_position"),
                       filters = gid.type,
                       values = gid,
                       mart = ensembl))
    g.list.human.hom.ext = as.data.table(getBM(attributes = c('ensembl_gene_id',
                                      'mmusculus_homolog_ensembl_gene',
                                      'mmusculus_homolog_associated_gene_name',
                                      'mmusculus_homolog_orthology_type',
                                      'mmusculus_homolog_orthology_confidence'),
                       filters = 'ensembl_gene_id',
                       values = g.list.human.ext$ensembl_gene_id,
                       mart = ensembl))
    #head(g.list.human.ext)
    #head(g.list.human.hom.ext)

    g.list.human.ext <- g.list.human.ext[g.list.human.hom.ext,,on='ensembl_gene_id']
    return(g.list.human.ext)
}


exp.human.file <- "baron-human.rds"
exp.mouse.file <- "baron-mouse.rds"
out.prefix <- "./baron"

sce.list <- list()
sce.list[["human"]] <- readRDS(exp.human.file)
sce.list[["mouse"]] <- readRDS(exp.mouse.file)

sce.list$human <- ssc.build(sce.list$human,display.name=rownames(sce.list$human))
sce.list$mouse <- ssc.build(sce.list$mouse,display.name=rownames(sce.list$mouse))
saveRDS(sce.list$human,sprintf("%s.sce.human.rds",out.prefix))
saveRDS(sce.list$mouse,sprintf("%s.sce.mouse.rds",out.prefix))

aid.names <- names(sce.list)
sce.list <- llply(names(sce.list),function(aid){
    sce <- sce.list[[aid]]
    sce$cell_type1 <- as.character(sce$cell_type1)
    ##hist(sce$total_features,breaks=50)
    expressed.ncell <- rowSums(assay(sce,"counts")>0)
    f.gene <- expressed.ncell > 2
    sce <- sce[f.gene,]
    return(sce)
})
names(sce.list) <- aid.names

sce.list$human$cell_type1[sce.list$human$cell_type1=="t_cell"] <- "T_cell"

table(sce.list$human$cell_type1)
table(sce.list$mouse$cell_type1)

#grep("^MT-",rownames(sce.list[["human"]]),value=T)
#grep("^Mt-",rownames(sce.list[["mouse"]]),value=T)

saveRDS(sce.list$human,sprintf("%s.sce.flt.human.rds",out.prefix))
saveRDS(sce.list$mouse,sprintf("%s.sce.flt.mouse.rds",out.prefix))

gid.mapping.df <- geneID.homolog(rownames(sce.list$human),gid.type="hgnc_symbol")

gid.mapping.slim <- unique(gid.mapping.df[,.(hgnc_symbol,
                                             mmusculus_homolog_associated_gene_name,
                                             mmusculus_homolog_orthology_type,
                                             mmusculus_homolog_orthology_confidence)])
gid.mapping.slim[,table(mmusculus_homolog_orthology_type, mmusculus_homolog_orthology_confidence)]
gid.mapping.slim <- gid.mapping.slim[mmusculus_homolog_orthology_type=="ortholog_one2one" & 
                                     mmusculus_homolog_orthology_confidence==1,]
gid.mapping.slim[,table(mmusculus_homolog_orthology_type, mmusculus_homolog_orthology_confidence)]

gid.mapping <- structure(gid.mapping.slim$mmusculus_homolog_associated_gene_name,
                         names=gid.mapping.slim$hgnc_symbol)

sce.list$human <- sce.list$human[names(gid.mapping),]
rownames(sce.list$human) <- gid.mapping[rownames(sce.list$human)]
rowData(sce.list$human)$feature_symbol <- rownames(sce.list$human)
rowData(sce.list$human)$display.name <- rownames(sce.list$human)

f.gene <- intersect(rownames(sce.list$human),rownames(sce.list$mouse))

sce.list$human <- sce.list$human[f.gene,]
sce.list$mouse <- sce.list$mouse[f.gene,]
rowData(sce.list$human)
rowData(sce.list$mouse)

saveRDS(sce.list$human,sprintf("%s.sce.flt.hom.human.rds",out.prefix))
saveRDS(sce.list$mouse,sprintf("%s.sce.flt.hom.mouse.rds",out.prefix))





