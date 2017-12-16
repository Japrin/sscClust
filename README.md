[![Travis CI Build Status](https://travis-ci.org/Japrin/sscClust.svg?branch=master)](https://travis-ci.org/Japrin/sscClust)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/bx9kcptlomn93auf/branch/master?svg=true)](https://ci.appveyor.com/project/Japrin/sscclust/branch/master)

# sscClust
simpler single cell RNAseq data clustering (sscClust), is a package implement mutilple functionalities which are basic procedures in single cell RNAseq data analysis, including variable genes identification, dimension reduction, clustering on reduced data. Also some new ideas were proposed, such as projecting data to a feature space using spearman correlation which make visualiaztion and clustering more efficient, clustering with subsampling and classfification which make it feasible to process thousands cells' data.

# Installation
Calculation will become expensive when the dataset become large, so we stongly recommand use R which linked to optimized BLAS library, such as ATLAS, MKL etc. For windows users, Microsoft R open is recommanded, and for Unix-alike users, please refer [this](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#BLAS) for how to compile R using external BLAS library.

To install this package, simply:
```
install.packages("devtools")
devtools::install_github("Japrin/sscClust")
```

# Example
Run the clustering pipeline for clustering using all data:
```
data("sce.Pollen")
sce.all <- ssc.run(sce.Pollen, subsampling=F, k.batch=5)
```
for clustering with subsampling:
```
data("sce.Pollen")
sce.sub <- ssc.run(sce.Pollen, subsampling=T, sub.frac = 0.8, k.batch=5)
```

More information can be found in the [vignette](http://htmlpreview.github.io/?https://github.com/Japrin/sscClust/blob/master/vignettes/sscClust.html).
