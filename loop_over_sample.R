rm(list = ls())

####Libraries#######
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(scds)
#devtools::install_github("keshavmot2/scanalysis")
library(scanalysis)
library(moments)
library(scran)
library(ggplot2)
library(PCAtools)
library(batchelor)
library(bluster)
library(pheatmap)
library(clustree)
library(Cairo)
library(BiocSingular)
library(cowplot)
library(batchelor)
library(dplyr)
library(GO.db)
library(topGO)
library(clusterProfiler)

setwd("/home/izem/Desktop/STAGE/")
bp.params <- MulticoreParam(workers = 5)

#sample.path_1_ctrl<-"/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_1_ctrl/outs/raw_feature_bc_matrix/"
sample.path_2_ctrl<-"/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_2_ctrl/outs/raw_feature_bc_matrix/"
sample.path_3_isch<-"/home/izem/Desktop/STAGE/Izem_Single-cell/Patient_IMI/IMI_data_1/outs/raw_feature_bc_matrix/"

sample.article_1<-"/home/izem/Desktop/STAGE/raw_feature_bc_matrix_1/"
sample.article_2<-"/home/izem/Desktop/STAGE/raw_feature_bc_matrix_2/"

#sce.sing_1_ctrl <- read10xCounts(sample.path_1_ctrl, col.names=TRUE, BPPARAM = bp.params)
sce.sing_2_ctrl<-read10xCounts(sample.path_2_ctrl, col.names=TRUE, BPPARAM = bp.params)
sce.sing_3_isch<-read10xCounts(sample.path_3_isch, col.names=TRUE, BPPARAM = bp.params)

#for the data of the article (Nature)
sce.sing_article_1<-read10xCounts(sample.article_1, col.names=TRUE, BPPARAM = bp.params)
sce.sing_article_2<-read10xCounts(sample.article_2, col.names=TRUE, BPPARAM = bp.params)

#for the data of the article (mouse)
sample.article_mice<-"/home/izem/Desktop/STAGE/raw_feature_bc_matrix_mouse/"
sce.sing_mice<-read10xCounts(sample.article_mice, col.names=TRUE, BPPARAM = bp.params)


all.sce<-list(sce.sing_2_ctrl,sce.sing_3_isch)

#all.sce<-list(sce.sing_article_1,sce.sing_article_2)



BiocManager::install("clusterProfiler")


#88888888888888888888888888888888888888888888888888888888888# EMPTY DROPS (start)
set.seed(1556)

i=0
empty_goutte<-list()
all.sce.emptied<-list()

for(n in all.sce){
  i=i+1
  current <-n
  empty_goutte[[i]]<-emptyDrops(counts(current),BPPARAM=bp.params)
  all.sce.emptied[[i]] <- current[,which(empty_goutte[[i]]$FDR <= 0.005)] 
  all.sce.emptied[[i]]<-cxds_bcds_hybrid(current)
  print(all.sce.emptied[[i]])
  
}

#88888888888888888888888888888888888888888888888888888888888# EMPTY DROPS (end)

#remove unexpressed genes

i=0
detected_genes<-list()
all.sce.detected<-list()
for(n in all.sce.emptied){
  i=i+1
  current <-n
  detected_genes[[i]]<-rowSums(counts(current)) > 0
  #all.sce.detected[[i]]<-current[detected_genes[[i]],]
  print(table(detected_genes[[i]]))
}


#888888888888888888888888888888888888888888888888888888888# Annotate genes (start)

i=0
genes<-list()
gene_annot<-list()
all.sce.annotate<-list()
ah <- AnnotationHub()
ens.mm.98 <- query(ah, c("Homo sapiens", "EnsDb", 98))[[1]] 
for(n in all.sce.detected){
  i=i+1
  current <-n
  genes[[i]] <- rowData(current)$ID
  gene_annot[[i]] <- AnnotationDbi::select(ens.mm.98, 
                                      keys = genes[[i]],
                                      keytype = "GENEID",
                                      columns = c("GENEID", "SEQNAME"))
  colnames(gene_annot[[i]])<-c("ID", "Chromosome")
  rowData(current) <- merge(rowData(current), gene_annot[[i]], by = "ID", sort=FALSE)
  rownames(rowData(current)) <- rowData(current)$ID
  all.sce.annotate[[i]]<-current
  print(rowData(current))
}


?AnnotationDbi::select

rowData(current)

#88888888888888888888888888888888888888888888888888888888# Annotate genes (end)


#8888888888888888888888888888888888888888888888888888888888888# Control Quality (start)


i=0
is.mito<-list()
all.sce.QCed<-list()
stats<-list()
for(n in all.sce.annotate){
  i=i+1
  current <-n
  is.mito[[i]] <- which(rowData(current)$Chromosome=="MT")
  stats[[i]]<-perCellQCMetrics(current, subsets=list(Mito=is.mito[[i]]), BPPARAM = bp.params)
  colData(current)<-cbind(colData(current),stats[[i]])
  current<-addPerFeatureQC(current,BPPARAM=bp.params)
  all.sce.QCed[[i]]<-current
}

#888888888888888888888888888888888888888888888888888888888888888# Control Quality (end)

all.sce.QCed

#888888888888888888888888888888888888888888888888888888888888888# Filtering cells (Start)

i=0
all.sce.filtered<-list()

for(n in all.sce.QCed){
  i=i+1
  current <-n
  low.lib <- isOutlier(current$sum,log = TRUE, type="lower")
  colData(current)$cell_sparsity <- 1 - (colData(current)$detected / nrow(current))
  rowData(current)$gene_sparsity <- (100 - rowData(current)$detected) / 100
  sparse.cells <- current$cell_sparsity > 0.990 #remove cells that express less than 0.005 genes
  #print(table(sparse.cells))
  #cat("table te3",i,sep="\n")
  mito.cells <- current$subsets_Mito_percent > 10
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


#table(sparse.cells)
#table(low.lib)

all.sce.filtered

#888888888888888888888888888888888888888888888888888888888888888# Filtering cells (End)

#888888888888888888888888888888888888888888888888888888888888888# Normalization + Dim Reduc (start)

i=0
set.seed(457)
all.sce.deconv<-list()
gene_var<-list()
hvgs<-list()
for(n in all.sce.filtered){
  i=i+1
  current <-n
  clust <- quickCluster(current, BPPARAM=bp.params)
  current <- computePooledFactors(current,
                                  clusters = clust,
                                  min.mean = 0.1,
                                  BPPARAM = bp.params)
  deconv.sf <- sizeFactors(current)
  current<-logNormCounts(current)
  gene_var[[i]]<- modelGeneVar(current)
  hvgs[[i]]<-getTopHVGs(gene_var[[i]], prop = 0.15)
  current<-scater::runPCA(current,subset_row=hvgs[[i]])
  current <- runUMAP(current,dimred="PCA") #stock UMAP
  current<-runTSNE(current,dimred="PCA") #stock tSNE
  rownames(current) <- uniquifyFeatureNames(rownames(current), rowData(current)$Symbol)
  all.sce.deconv[[i]]<-current
}

#all.sce<-lapply(all.sce, logNormCounts)
#all.dec <- lapply(all.sce.deconv,modelGeneVar)
#all.hvgs <- lapply(all.dec, getTopHVGs, prop=0.15)


current <- runUMAP(current)
plotReducedDim(current, dimred = "UMAP", colour_by = "detected")


pca_pct_variance <- data.frame(variance = attr(reducedDim(all.sce.deconv[[2]], "PCA"), "percentVar"))
variance<-pca_pct_variance$variance

plotReducedDim(all.sce.deconv[[2]], dimred = "PCA", colour_by = "detected")
plotReducedDim(all.sce.deconv[[2]],dimred = "UMAP", colour_by = "label")

current<-all.sce.deconv[[2]]

current$label<-clust


#888888888888888888888888888888888888888888888888888888888888888# Normalization + Dim Reduc (end)


#88# Data integration (Start) (BROAD)#########

#i've already done it but i create the all.dec variable in order to continue with the variables of the workshop
all.dec <- lapply(all.sce.deconv, modelGeneVar)

universe <- Reduce(intersect, lapply(all.sce.deconv, rownames))

# Subsetting the SingleCellExperiment object.
uni.sce <- lapply(all.sce.deconv, function(x){x[universe,]})
# Also subsetting the variance modelling results, for convenience.
uni.dec <- lapply(all.dec, function(x){x[universe,]})


# Renormalizing to adjust for differences in depth.
normed.sce <-rescaled <- do.call(multiBatchNorm, uni.sce)




# Identifying a set of HVGs using stats from all batches.
combined.dec <- combineVar(uni.dec)
combined.hvg <- getTopHVGs(combined.dec, n=5000)
#I can choose a big set of HVGs (for ex 5000) or I can choose all the genes with the biological component > 0 
chosen.hvgs <- combined.dec$bio > 0
combined.dec$chosenHvg <- chosen.hvgs




##888888# linear regression 

linear_rescaled <- rescaleBatches(rescaled)

##888888# linear regression 

### It's really just for the code --------
data_1_ctrl<-normed.sce[[1]]
data_2_ctrl<-normed.sce[[2]]

rowData(data_1_ctrl) <- rowData(data_2_ctrl)

data_1_ctrl$batch<-"ctrl_1"
data_2_ctrl$batch<-"ctrl_2"

uncorrected <- cbind(data_1_ctrl,data_2_ctrl)

####999999999999999##

set.seed(1000101002)
mnn.out <- do.call(fastMNN, c(rescaled, 
                                  list(subset.row=chosen.hvgs, BSPARAM=RandomParam())))

#save it in a variable because he desapear with the function SNNgraph
batch<-mnn.out$batch

metadata(mnn.out)$merge.info$lost.var

mnn.out.corre.dim <- dim(reducedDim(mnn.out, "corrected"))
mnn.out.recon.dim <- dim(assay(mnn.out, "reconstructed"))

#g <- buildSNNGraph(mnn.out, use.dimred="corrected")
#colLabels(mnn.out) <- factor(igraph::cluster_louvain(g)$membership)


# add feature selection outcome to mmn.out
# used in other analyses later.
columnsToGet <- setdiff(colnames(combined.dec), "per.block")
combined.dec.df <- combined.dec[,columnsToGet] %>%
  data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  filter(ID %in% rownames(rowData(mnn.out)))
rotationMat <- rowData(mnn.out)$rotation
rowData(mnn.out)$ID <- rownames(rowData(mnn.out))
rowData(mnn.out)$rotation <- NULL


rowData(mnn.out) <- rowData(mnn.out) %>%
  data.frame() %>%
  left_join(combined.dec.df, by="ID") %>%
  DataFrame()


# add rotation back
rowData(mnn.out)$rotation <- rotationMat

#in order to skip the uncorrected
rowData(mnn.out)$Symbol<-rownames(mnn.out)



#tidy
rm(columnsToGet, rotationMat) 

colDataList <- lapply(rescaled, function(x){colData(x)})
colDataDf <- do.call(rbind, colDataList)
colData(mnn.out) <- DataFrame(colDataDf)

snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected", k=20)
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
colLabels(mnn.out) <- factor(clusters.mnn)
clusterTab.mnn <- data.frame(clusters=clusters.mnn, batch=batch, source=mnn.out$label)
#clusters and source col are the same 


ClusterInfo.mnn <- clusterTab.mnn %>% 
  as.data.frame() %>%
  group_by(clusters,batch) %>%
  summarise(cells = n())


#to see the number of cells involved in each row
mp1 <- ggplot(data=ClusterInfo.mnn, aes(x=clusters,y=cells, fill=batch)) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  geom_col() +
  theme(legend.text = element_text(size = 7))

plot_grid(mp1,ncol=1)

mnn.out$batch<-batch

set.seed(10101010)
mnn.out <- runTSNE(mnn.out, dimred="corrected")
gridExtra::grid.arrange(
  plotTSNE(mnn.out, colour_by="label", text_by="label", text_colour="red"),
  plotTSNE(mnn.out, colour_by="batch")
)

#the % of variance lost 
round(metadata(mnn.out)$merge.info$lost.var,2)


## T-SNE with interesting genes (just for the code) -----------

genesToShow <- c("PAX7","PECAM1","CDKN2A","APOE")


tmpInd <- which(rownames(mnn.out) %in% genesToShow)


p <- plotTSNE(mnn.out, colour_by = genesToShow[1], by_exprs_values="reconstructed",shape_by = "batch")
p <- p + ggtitle(paste("Cellules SAT", genesToShow[1]))
psat <- p

p <- plotTSNE(mnn.out, colour_by = genesToShow[2], by_exprs_values="reconstructed",shape_by = "batch")
p <- p + ggtitle(paste("Cellules Endothéliales", genesToShow[2]))
pendo <- p


p <- plotTSNE(mnn.out, colour_by = genesToShow[3], by_exprs_values="reconstructed",shape_by = "batch")
p <- p + ggtitle(paste("Super Cellules SAT", genesToShow[3]))
pSsat <- p

p <- plotTSNE(mnn.out, colour_by = genesToShow[4], by_exprs_values="reconstructed",shape_by = "batch")
p <- p + ggtitle(paste("Monocytes", genesToShow[4]))
proli <- p


plot(proli)
gridExtra::grid.arrange(psat + theme(legend.position="left"),
                        pendo + theme(legend.position="left"),
                        pSsat + theme(legend.position = "left"),
                        ncol=2)


sum(all.sce.QCed[[2]]$sum)
## Comparing with the uncorrected one (just for the code) ----------------

set.seed(1111001)
uncorrected <- runTSNE(uncorrected, dimred="PCA")
plotTSNE(uncorrected, colour_by="batch")

p <- plotTSNE(uncorrected, colour_by = genesToShow[1])
p <- p + ggtitle(paste("Cellules SAT", genesToShow[1]))
psat <- p

p <- plotTSNE(uncorrected, colour_by = genesToShow[2])
p <- p + ggtitle(paste("Cellules Endothéliales", genesToShow[2]))
pendo <- p

p <- plotTSNE(uncorrected, colour_by = genesToShow[3])
p <- p + ggtitle(paste("Super Cellules SAT", genesToShow[3]))
pSsat <- p

p <- plotTSNE(uncorrected, colour_by = genesToShow[4])
p <- p + ggtitle(paste("Cellules en prolifération", genesToShow[4]))
proli <- p

gridExtra::grid.arrange(psat + theme(legend.position="bottom"),
                        pendo + theme(legend.position="bottom"),
                        pSsat + theme(legend.position = "bottom"),
                        proli + theme(legend.position = "bottom"),
                        ncol=2)

#999999999999999999999999999999# Batch correction Diagnosis

ariVec <- vector(mode = "numeric", length = 2)
sampleNames <- c("isch","ctrl")
names(ariVec) <- c("isch","ctrl")

for (i in 1:2) {
  ariVec[i] <- pairwiseRand(
    ref=as.integer(colLabels(rescaled[[i]])),
    alt=as.integer(clusters.mnn[linear_rescaled$batch==sampleNames[i]]),
    mode="index")
}
ariVec <- round(ariVec,2)
ariVec




#888888888888888888888888888888888888888888888888888888888888888# Data integration (end)


#888888888888888888888888888888888888888888888888888888888888888# Differential expression analysis (start)




#888888888888888888888888888888888888888888888888888888888888888# Differential expression analysis (end)



