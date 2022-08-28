import os

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects.numpy2ri as rpyn
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri as pdri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import StrVector
import scanpy as sc
import scipy , pandas,numpy as np

'''###This script has one goal , suppose someone wants to begin with raw matrix counts , i can't begin with scanpy so i have to begin with R and DropletUtils ect...

importr('BiocParallel')
importr('DropletUtils')
importr('scater')
#variable=robjects.r('sample.path_1 <- "/home/izem/Desktop/STAGE/Izem_Single-cell/Patient_IMI/IMI_data_1/outs/raw_feature_bc_matrix/"')

robjects.r('sample.path_1 <- "/home/izem/Desktop/STAGE/Izem_Single-cell/Patient_IMI/IMI_data_1/outs/raw_feature_bc_matrix/"')
#print(robjects.r['sample.path_1'])
robjects.r('bp.params <- MulticoreParam(workers = 5)')

dataa=robjects.r('sce.sing <- read10xCounts(sample.path_1, col.names=TRUE, BPPARAM = bp.params)')

dolto=robjects.r('set.seed(1556);e.out <- emptyDrops(counts(sce.sing),BPPARAM=bp.params);sce.pbmc <- sce.sing[,which(e.out$FDR <= 0.005)] ;sce.pbmc<-logNormCounts(sce.pbmc)')
#il garde que le dernier


robjects.r('explore<-t(as.matrix(logcounts(sce.pbmc)));very_large<-list();for(x in 1:length(colnames(sce.pbmc))){xxx<-unname(explore[x,]);very_large[[x]]<-xxx};large<-unlist(very_large)')
#what i do in this R script:
# create a transposed matrix (because annData have feature on cols and cell on rows) 
# the rest is combining all the values of this matrix in a large vector

#i'm force to first convert the number into numpy and then convert it to int
dim_row=int(rpyn.rpy2py(robjects.r('dim(explore)[1]')))
dim_col=int(rpyn.rpy2py(robjects.r('dim(explore)[2]')))


#minus=np.asarray(robjects.r['large'])

minus=rpyn.rpy2py(robjects.r['large'])

#print(type(int(rpyn.rpy2py(dim_row))))


cbon_douka=minus.reshape((dim_row,dim_col))

#var_a=robjects.r('exprs <- t(as.matrix(logcounts(sce.pbmc)))') #i have to transform it to a matrix i think
var_b=pdri.rpy2py(robjects.r('col_data <- as.data.frame(colData(sce.pbmc))'))
var_c=pdri.rpy2py(robjects.r('row_data <- as.data.frame(rowData(sce.pbmc))'))

#print(type(var_b))

#print(robjects.r['sce.pbmc'])

adata_sce = sc.AnnData(X = cbon_douka, obs = var_b, var = var_c)
'''


'''importr('DropletUtils') ; importr('scran') ; importr("scuttle") ;importr("SingleCellExperiment") ; importr("scater") ; importr("anndata")


robjects.r("to_sce<-anndata::read_h5ad('/home/izem/test.h5ad')")

robjects.r("sce_janin<-SingleCellExperiment(assays=list(counts=t(to_sce$X)),rowData=list(to_sce$var),colData=list(to_sce$obs))")

robjects.r("libsize=unname(librarySizeFactors(sce_janin))")

robjects.r("set.seed(457);clust <- quickCluster(sce_janin);sce_janin <- computePooledFactors(sce_janin, clusters = clust, min.mean = 0.1)")

robjects.r("deconv.sf <- sizeFactors(sce_janin);sce_janin <- logNormCounts(sce_janin)")

robjects.r("colnames(colData(sce_janin))[7]<-'deconvSizeFactor';sce_janin$libSizeFactor=libsize")

robjects.r("from_sce=anndata::AnnData(X=t(logcounts(sce_janin)),obs=as.data.frame(colData(sce_janin)),var =as.data.frame(rowData(sce_janin)))")'''

#script="anndata::write_h5ad(from_sce,"+os.getcwd()+"/complete.h5ad"

#robjects.r(script)

import MNN

MNN.convertR("/home/izem/test.h5ad")
