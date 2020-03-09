.onLoad <- function(libname , pkgname) {
	if(getRversion() >= "2.15.1"){
		utils::globalVariables(c(".",".N",".SD","variable","value","Group","N","x","y",
								 "geneID","aid","cellID","cellID.raw","dataset",
								 "ClusterID", "F.adjp", "F.pvalue", "HSD.padj.min", "HSD.padj.min.diff",
								 "P.Value","adj.P.Val","logFC",
								 "o.Exp","bin.Exp",
								 "dprime","vardprime","freq.sig",
								 "silhouette","silhouette.rank",
								 "median.F.rank","median.rank","meta.cluster", "meta.cluster.RF", "meta.cluster.raw",
								 "RF.meta.cluster.max", "RF.meta.cluster.sec", "RF.prob.max RF.prob.sec",
								 "RF.class.max", "RF.class.sec", "RF.prob.max", "RF.prob.sec",
								 "pred","actual"
								 ))
	}
}
