

library(MouseGastrulationData)
sce.chimera <- WTChimeraData(samples=5:10)
sce.chimera
rownames(sce.chimera) <- uniquifyFeatureNames(
  rowData(sce.chimera)$ENSEMBL, rowData(sce.chimera)$SYMBOL)

drop <- sce.chimera$celltype.mapped %in% c("stripped", "Doublet")
sce.chimera <- sce.chimera[,!drop]

sce.chimera <- logNormCounts(sce.chimera)

dec.chimera <- modelGeneVar(sce.chimera, block=sce.chimera$sample)
chosen.hvgs <- dec.chimera$bio > 0

par(mfrow=c(1,2))
blocked.stats <- dec.chimera$per.block
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(current$mean, current$total, main=i, pch=16, cex=0.5,
       xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit <- metadata(current)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
}


library(batchelor)
set.seed(01001001)
merged <- correctExperiments(sce.chimera, 
                             batch=sce.chimera$sample, 
                             subset.row=chosen.hvgs,
                             PARAM=FastMnnParam(
                               merge.order=list(
                                 list(1,3,5), # WT (3 replicates)
                                 list(2,4,6)  # td-Tomato (3 replicates)
                               )
                             )
)


metadata(merged)$merge.info$lost.var

g <- buildSNNGraph(merged, use.dimred="corrected")
clusters <- igraph::cluster_louvain(g)
colLabels(merged) <- factor(clusters$membership)

merged <- runTSNE(merged, dimred="corrected", external_neighbors=TRUE)
merged <- runUMAP(merged, dimred="corrected", external_neighbors=TRUE)


# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.
summed <- aggregateAcrossCells(merged, 
                               id=colData(merged)[,c("celltype.mapped", "sample")])

summed$ncells

label <- "Mesenchyme"
current <- summed[,label==summed$celltype.mapped]


# Creating up a DGEList object for use in edgeR:
library(edgeR)
y <- DGEList(counts(current)[,1:4], samples=colData(current)[1:4,])
#y <- DGEList(counts(current), samples=colData(current))
y

discarded <- current$ncells < 10
y <- y[,!discarded]
summary(discarded)


keep <- filterByExpr(y, group=current$tomato)
y <- y[keep,]
summary(keep)

y <- calcNormFactors(y)
y$samples

design <- model.matrix(~factor(tomato), y$samples)
design

y <- estimateDisp(y, design)
summary(y$trended.dispersion)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)

res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))

topTags(res)

?glmQLFTest
#####now do it with findMarkers######

big_index<-c(which(merged$batch=="5" & merged$sample=="5" & merged$tomato==TRUE & merged$pool=="3" & merged$celltype.mapped=="Mesenchyme")
             ,which(merged$batch=="6" & merged$sample=="6" & merged$tomato==FALSE & merged$pool=="3" & merged$celltype.mapped=="Mesenchyme")
             ,which(merged$batch=="7" & merged$sample=="7" & merged$tomato==TRUE & merged$pool=="4" & merged$celltype.mapped=="Mesenchyme")
             ,which(merged$batch=="8" & merged$sample=="8" & merged$tomato==FALSE & merged$pool=="4" & merged$celltype.mapped=="Mesenchyme"))

nahi<-merged[,big_index]

#now let's validate our technique

markers_default <- findMarkers(
  nahi,
  assay.type="logcounts",
  groups = factor(nahi$tomato), # clusters to compare
  test.type = "t",   # t-test (default)
  direction = "any", # test for either higher or lower expression (default)
  lfc = 0, # null hypothesis log-fold-change = 0 (default)
  pval.type = "any" # ranking of p-values based on any comparison (default)
)


View(cbind(as.data.frame(markers_default[[2]]),rownames(markers_default[[2]])))

#i found almost the same DE genes . The topTags(res) ie top 10 were found in the top 30 of the marker output 
#(except for the last one Prps2 which was at 101 position) 
#I also tested with Allantois cells and I had more than 5 of topTags(res) found in top 20 or 30






