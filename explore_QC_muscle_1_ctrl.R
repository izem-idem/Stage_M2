rm(list = ls())

library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(scds)
library(Seurat)
library(SeuratObject)
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
  
sample.path_1 <- "/home/izem/Desktop/STAGE/Izem_Single-cell/Patient_IMI/IMI_data_1/outs/raw_feature_bc_matrix/"

sample.path_1_ctrl<-"/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_1_ctrl/outs/raw_feature_bc_matrix/"


sce.sing <- read10xCounts(sample.path_1_ctrl, col.names=TRUE, BPPARAM = bp.params)


#---------------------
#to try something on the IM one
sce.sing<-sce.sing_3_isch

set.seed(1556)
e.out <- emptyDrops(counts(sce.sing),BPPARAM=bp.params)
sce.pbmc <- sce.sing[,which(e.out$FDR <= 0.005)] #on récupére les indices des barcodes qui ont un FDR inférieur ou égal à 0.005
unfiltered <- sce.pbmc
sce.pbmc



#--------------------------------

detected_genes <- rowSums(counts(sce.pbmc)) > 0

table(detected_genes)

sce.pbmc <- sce.pbmc[detected_genes,]

#-------------------------------------------


ah <- AnnotationHub()
ens.mm.98 <- query(ah, c("Homo sapiens", "EnsDb", 98))[[1]] 

genes <- rowData(sce.pbmc)$ID
gene_annot <- AnnotationDbi::select(ens.mm.98, 
                                    keys = genes,
                                    keytype = "GENEID",
                                    columns = c("GENEID", "SEQNAME"))

colnames(gene_annot)<-c("ID", "Chromosome")
rowData(sce.pbmc) <- merge(rowData(sce.pbmc), gene_annot, by = "ID", sort=FALSE)
rownames(rowData(sce.pbmc)) <- rowData(sce.pbmc)$ID


#----------------------------




is.mito <- which(rowData(sce.pbmc)$Chromosome=="MT")
sce.pbmc<- addPerCellQC(sce.pbmc, subsets=list(Mito=is.mito), BPPARAM = bp.params)


#View(as.data.frame(colData(sce.final)[3:4]))



#--------------------------------------------------



#hist(sce.pbmc$detected[which(sce.pbmc$detected<200)],breaks = 100)



#colData(sce.pbmc)$log10_counts<-log10(sce.pbmc$sum+1)

low.lib <- isOutlier(sce.pbmc$sum, type="lower",log = TRUE)


isOutlier(sce.pbmc$detected,log = TRUE ,type = "higher")

table(low.lib)

#------------------------------------------

sce.pbmc <- addPerFeatureQC(sce.pbmc, BPPARAM = bp.params)
rowData(sce.pbmc)

colData(sce.pbmc)$cell_sparsity <- 1 - (colData(sce.pbmc)$detected / nrow(sce.pbmc))
rowData(sce.pbmc)$gene_sparsity <- (100 - rowData(sce.pbmc)$detected) / 100

sparse.cells <- sce.pbmc$cell_sparsity > 0.990
mito.cells <- sce.pbmc$subsets_Mito_percent > 10
min.cells <- 1 - (20 / ncol(sce.pbmc))
sparse.genes <- rowData(sce.pbmc)$gene_sparsity > min.cells

table(sparse.genes)
table(mito.cells)

table(low.lib,sparse.cells)

table(low.lib,mito.cells)

table(sparse.cells)

#sce.pbmc <- cxds_bcds_hybrid(sce.pbmc)
#pour detecter les doublets

#----------------------Removing low Q cells


da_one<-which(sparse.cells==FALSE)
da_two<-which(mito.cells==FALSE)
da_three<-which(low.lib==FALSE)

inter_one<-intersect(da_one,da_two)
final_inter<-intersect(inter_one,da_three)

str(final_inter)

#645 cells

sce.final<-sce.pbmc[,final_inter]

sce.final

##########NORMALISATION (Seurat) -cancelled################

#from sce to seurat in order to do the normalisation
seu.resume<-sce_to_seurat(sce.final)

Seurat::GetAssayData(seu.resume,slot="counts")

seu.resume <- Seurat::NormalizeData(seu.resume,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)


##test classic workflow
seu.resume <- FindVariableFeatures(seu.resume, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seu.resume)
seu.resume <- ScaleData(seu.resume, features = all.genes)

seu.resume <- RunPCA(seu.resume, features = VariableFeatures(object = seu.resume))
###



#stored in seu.resume[["RNA"]]@data 

summary(seu.resume[["RNA"]]@data)

hist(summary(seu.resume[["RNA"]]@data)[,3])

###Explore normalization for the 10 firsts cells##
dfeyme<-as.data.frame(summary(seu.resume[["RNA"]]@data)[,-1])

dfeyme<-dfeyme[seq(1,26206,1),]

boxplot(dfeyme[,2]~dfeyme[,1],col=seq(1,15,1))
####


###Deconvolution##############

set.seed(457) # clusters with PCA from irlba with approximation
clust <- quickCluster(sce.final, BPPARAM=bp.params) # slow with all cells.
table(clust)

sce.final <- computePooledFactors(sce.final,
                            clusters = clust,
                            min.mean = 0.1,
                            BPPARAM = bp.params)
deconv.sf <- sizeFactors(sce.final)
summary(deconv.sf)


sce.final <- logNormCounts(sce.final)


print(assayNames(sce.final))
logcounts(sce.final)

tableau<-as.data.frame(summary(logcounts(sce.final))[seq(1,26206,1),-1])
#check the 15 first cells 

boxplot(tableau[,2]~tableau[,1],col=seq(1,15,1))

par(mfrow=c(1,2))


######SCT transform (SEURAT) -cancelled###########


#renames cols so that the seurat's functions works
colnames(seu.resume@meta.data)[2] <-"nCount_RNA"

colnames(seu.resume@meta.data)[3] <-"nFeature_RNA"


seu.resume <- PercentageFeatureSet(seu.resume, pattern = "^MT-", col.name = "percent.mt")

seu.resume <- SCTransform(seu.resume, verbose = TRUE)

tabla<-as.data.frame(summary(seu.resume@assays[["SCT"]]@counts)[seq(1,23779,1),-1])

boxplot(tabla[,2]~tabla[,1],col=seq(1,15,1))

############DIMENSIONNALITY REDUCTION (SEURAT) - cancelled ##########


seu.resume <- Seurat::RunPCA(seu.resume)

Seurat::DimPlot(seu.resume, reduction = "pca")

##To see the variance explained by each PCs : 

mat <- Seurat::GetAssayData(seu.resume, assay = "RNA", slot = "scale.data")
pca <- seu.resume@reductions[["pca"]]

total_variance <- sum(matrixStats::rowVars(mat))
varExplained = eigValues / total_variance





#####COntinue with deconvol#####

rownames(sce.final) <- uniquifyFeatureNames(rownames(sce.final), rowData(sce.final)$Symbol)

gene_var <- modelGeneVar(sce.final)

# plot the variance against mean of expression with fitted trend line
gene_var %>% 
  # convert to tibble for ggplot
  as_tibble() %>% 
  # make the plot
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_line(aes(y = tech), colour = "dodgerblue", size = 1) +
  labs(x = "Mean of log-expression", y = "Variance of log-expression")


hvgs<-getTopHVGs(gene_var, prop = 0.15)

sce.final <- runPCA(sce.final, subset_row = hvgs)

sce.final<-runTSNE(sce.final,dimred="PCA")

g <- buildSNNGraph(sce.final, k=5, use.dimred='PCA') # pour k=5, 10c / k=4 11c / k=3, 14c / k=2 19c 
clustrr <- igraph::cluster_louvain(g)$membership
sce.final$label  <- factor(clustrr)



pca_pct_variance <- data.frame(variance = attr(reducedDim(sce.final, "PCA"), "percentVar"))
pca_pct_variance$PC <- 1:nrow(pca_pct_variance)


table(clustrr)

# visualise percentage variance explained by PCs (scree plot)
pca_pct_variance %>% 
  ggplot(aes(PC, variance)) +
  geom_col() 
  labs(y = "Variance explained (%)")

plotReducedDim(sce.final, dimred = "TSNE", colour_by = "CKM",text_by = "label")

#subset the sce.final to keep only cluster 7 cells

sce.sat.final<-sce.final[,which(sce.final$label=="7")]

plotReducedDim(sce.sat.final,dimred = "TSNE",colour_by="CDKN2A")


#plotReducedDim(sce.final, dimred = "PCA", ncomponents = 10, colour_by = "detected")

ggcells(sce.final, aes(x = PCA.1, y = PCA.2)) +
  geom_point(size = 0.5) +
  labs(x = "PC1", y = "PC2", colour = "Sample")


###compare the sat cells ctrl vs sat cells isch####

ctrl.sce<-all.sce.deconv[[2]]


table(clust)

ctrl.sce<-ctrl.sce[,which(ctrl.sce$label=="10")]

plotReducedDim(ctrl.sce, dimred = "TSNE", colour_by = "CDKN2A",text_by = "label")+ggtitle("Patient_ctrl")

#In order to see a difference.
gridExtra::grid.arrange(
  plotReducedDim(ctrl.sce, dimred = "TSNE", colour_by = "CDKN2A")+ggtitle("Patient_ctrl"),
  plotReducedDim(sce.sat.final,dimred = "TSNE",colour_by="CDKN2A")+ggtitle("Patient_isch")
)

##which marker for which cluster ! #######

marker.isch <- scoreMarkers(sce.final,sce.final$label)

chosen <- marker.isch[["11"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
#sdd<-as.data.frame(cbind(ordered[,1:4],rownames(ordered))) # showing basic stats only, for brevity.
ordered[1:10,1:4]

#to test the marker
plotReducedDim(sce.final, dimred = "TSNE",text_by = "label")+geom_point(aes(color=factor(sce.final$label)))


table(clust)

#which position
which(rownames(ordered)=="CD52")

#let's try SingleR#####

library(celldex)
ref.data <- HumanPrimaryCellAtlasData(ensembl=FALSE)

library(SingleR)
predictions <- SingleR(test=sce.final, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main,clusters =sce.final$label)

chafou<-as.data.frame(predictions)



#####SCT transform (BROAD) -cancelled############

set.seed(44) # for reproducibility
vst_out <- sctransform::vst(
  counts(sce.final), # A matrix of UMI counts with genes as rows and cells as columns
  latent_var = c('log_umi'), # The independent variables to regress out as a character vector
  return_gene_attr = TRUE, # Make cell attributes part of the output
  return_cell_attr = TRUE, # Calculate gene attributes and make part of output
  verbosity = 0 # An integer specifying what to show (0: nothing, 1: messages, 2: + progress bar)
)


vssctransform::plot_model_pars(
  vst_out, # The output of a vst run
  verbosity = 1 # Messages only, no progress bar
)


geneOverlap <- rownames(sce.final) %in% rownames(vst_out$y)
if(!all(geneOverlap))
{
  table(rownames(sce.final) %in% rownames(vst_out$y))
  tmpInd <- which(rownames(sce.final) %in% rownames(vst_out$y))
  sce.final <- sce.final[tmpInd,]
  assayNames(sce.final)
}

vstMat <- as(vst_out$y[rownames(sce.final),], "dgCMatrix")
all(colnames(vstMat) == sce.final$Barcode)
all(rownames(vstMat) == rownames(sce.final))

assay(sce.final, "sctrans_norm", withDimnames=FALSE) <- vstMat

#Feature selection 

gene_var <- modelGeneVar(sce.final)

vary_genes<-getTopHVGs(gene_var, n=1000)

sce.final <- scater::runPCA(sce.final, subset_row = vary_genes)
pca_pct_variance <- data.frame(variance = attr(reducedDim(sce.final, "PCA"), "percentVar"))

###PCA DIAGNOSTIC (BROAD)

explan_pcs <- getExplanatoryPCs(sce.final,
                                variables = c(
                                  "sum",
                                  "detected",
                                  "subsets_Mito_percent"
                                )
)


plotExplanatoryPCs(explan_pcs/100)

plotExplanatoryVariables(sce.final,
                         variables = c(
                           "sum",
                           "detected",
                           "subsets_Mito_percent"
                         ))

chosen_elbow <- findElbowPoint(pca_pct_variance$variance)
chosen_elbow


subset.row=NULL
sce.final <- denoisePCA(sce.final, technical = gene_var,subset.row)
ncol(reducedDim(sce.final, "PCA"))


sum(variance[1:5])
