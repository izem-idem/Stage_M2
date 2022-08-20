
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(scds)
library(scanalysis)
library(scran)
library(ggplot2)
library(dplyr)
library(scDblFinder)
library(stringr)

bp.params <- MulticoreParam(workers = 5)

var.path.a<-"/home/izem/Desktop/Data_2/MTX_cKO_Nonane/outs/raw_feature_bc_matrix/"

var.path.b<-"/home/izem/Desktop/Data_2/MTX_cKO_TCDD/outs/raw_feature_bc_matrix/"

var.path.c<-"/home/izem/Desktop/Data_2/MTX_CTRL_Nonane/outs/raw_feature_bc_matrix/"

var.path.d<-"/home/izem/Desktop/Data_2/MTX_CTRL_TCDD/outs/raw_feature_bc_matrix/"


sce.nuc.cko_NONANE<-read10xCounts(var.path.a)
sce.nuc.cko_TCDD<-read10xCounts(var.path.b)
sce.nuc.ctrl_NONANE<-read10xCounts(var.path.c)
sce.nuc.ctrl_TCDD<-read10xCounts(var.path.d)


all.sce<-list(sce.nuc.cko_NONANE,sce.nuc.cko_TCDD,sce.nuc.ctrl_NONANE,sce.nuc.ctrl_TCDD)

set.seed(1556)

i=0
empty_goutte<-list()
all.sce.emptied<-list()

for(n in all.sce){
  i=i+1
  current <-n
  empty_goutte[[i]]<-emptyDrops(counts(current),BPPARAM=bp.params)
  all.sce.emptied[[i]] <- current[,which(empty_goutte[[i]]$FDR <= 0.005)] 
  #all.sce.emptied[[i]]<-cxds_bcds_hybrid(current) #for the article dataset 
  print(all.sce.emptied[[i]])
  
}

#######################################################################################

i=0
detected_genes<-list()
all.sce.detected<-list()
for(n in all.sce.emptied){
  i=i+1
  current <-n
  detected_genes[[i]]<-rowSums(counts(current)) > 0
  all.sce.detected[[i]]<-current[detected_genes[[i]],]
  print(table(detected_genes[[i]]))
}

####################################################################################

i=0
is.mito<-list()
all.sce.QCed<-list()
stats<-list()
for(n in all.sce.detected){
  i=i+1
  current <-n
  is.mito[[i]] <- which(startsWith(rowData(current)$Symbol,"MT"))
  stats[[i]]<-perCellQCMetrics(current, subsets=list(Mito=is.mito[[i]]), BPPARAM = bp.params)
  colData(current)<-cbind(colData(current),stats[[i]])
  current<-addPerFeatureQC(current,BPPARAM=bp.params)
  all.sce.QCed[[i]]<-current
}


###################################################################################



i=0
all.sce.filtered<-list()

for(n in all.sce.QCed){
  i=i+1
  current <-n
  low.lib <- isOutlier(current$sum,log = TRUE, type="lower")
  print(attr(low.lib, "thresholds")[1])
  colData(current)$cell_sparsity <- 1 - (colData(current)$detected / nrow(current))
  rowData(current)$gene_sparsity <- (100 - rowData(current)$detected) / 100
  sparse.cells <- current$cell_sparsity > 0.990 #remove cells that express less than 0.01 genes ie 1 % 
  #print(table(sparse.cells))
  #cat("table te3",i,sep="\n")
  mito.cells <- current$subsets_Mito_percent > 0.1
  min.cells <- 1 - (20 / ncol(current)) 
  sparse.genes <- rowData(current)$gene_sparsity > min.cells # remove genes that are expressed in less than 20 cells
  da_one<-which(sparse.cells==FALSE)
  da_two<-which(mito.cells==FALSE)
  da_three<-which(low.lib==FALSE) # vu le nombre faible de cellules je ne vais pas utiliser ce critére pour purifier
  inter_one<-intersect(da_one,da_two)
  final_inter<-intersect(inter_one,da_three)
  current<-current[,final_inter]
  all.sce.filtered[[i]]<-current
}

all.sce.filtered

# ok jusque là , besah pour le doublet detection je vais devoir choisir si je veux le faire ou pas 

#cxds_bcds_hybrid

all.sce.filtered[[1]]<-cxds_bcds_hybrid(all.sce.filtered[[1]])
all.sce.filtered[[2]]<-cxds_bcds_hybrid(all.sce.filtered[[2]])
all.sce.filtered[[3]]<-cxds_bcds_hybrid(all.sce.filtered[[3]])
all.sce.filtered[[4]]<-cxds_bcds_hybrid(all.sce.filtered[[4]])

table(all.sce.filtered[[2]]$hybrid_score>1)


