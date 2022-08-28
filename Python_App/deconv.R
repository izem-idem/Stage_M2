library(anndata)
library(scater)
library(scran)
library(SingleCellExperiment)
library(DropletUtils)
library(scuttle)

args <- commandArgs(TRUE)

#remotes::install_github("rstudio/reticulate")

setwd(args[1])

to_sce<-anndata::read_h5ad(args[2])


#to_sce<-anndata::read_h5ad("/home/izem/PycharmProjects/Explore/test.h5ad")



#save[450:460,3000:3010] #en comparant had le subset avec l'anndata de scanpy c'est la même chose du coup on est bon 

#TO SCE REALLY NOW

sce_janin<-SingleCellExperiment(assays=list(counts=t(to_sce$X)),rowData=list(to_sce$var),colData=list(to_sce$obs))

sce_janin

libsize=unname(librarySizeFactors(sce_janin))

set.seed(457) # clusters with PCA from irlba with approximation
clust <- quickCluster(sce_janin) # slow with all cells.
table(clust)

sce_janin <- computePooledFactors(sce_janin,
                                  clusters = clust,
                                  min.mean = 0.1)
                                  
deconv.sf <- sizeFactors(sce_janin)
summary(deconv.sf)

hist(sce_janin$sizeFactor)

sce_janin <- logNormCounts(sce_janin)

#plot(libsize, deconv.sf, xlab="Library size factor",ylab="Deconvolution size factor", log='xy', pch=16)

sce_janin$deconv=deconv.sf
sce_janin$libSizeFactor=libsize



#back to anndata NOW ! 

from_sce=anndata::AnnData(X=t(logcounts(sce_janin)),obs=as.data.frame(colData(sce_janin)),var =as.data.frame(rowData(sce_janin)))

#from_sce$X[450:460,3000:3010]

anndata::write_h5ad(from_sce,"complete.h5ad")




