#First we'll create the famous uncorrected SCE object.

#you have to first run the script loop_over_sample to create the variable all.sce.annotate 
#run it up to all.sce.deconv
#we start from all.sce.deconv and we do the following: 
i=0
for (n in all.sce.deconv) {
  i=i+1
  current<-n
  g <- buildSNNGraph(n, k=10, use.dimred='PCA')
  clust <- igraph::cluster_louvain(g)$membership
  current$label  <- factor(clust)
  all.sce.deconv[[i]]<-current
}

all.sce.deconv[[2]]$leiden




#le groupe 1 le cluster mlih c'est louvain , le groupe 2 c'est mieux de faire le cluster leiden
#mais louvain il est bien chwiya pour le groupe 2 car il rate que le cluster 9 (myocytes)

table(all.sce.deconv[[2]]$label)

#even though we already have the information regarding the model gene var in the all.sce.deconv loop
#we repeat it to cary on with the workshop
all.dec <- lapply(all.sce.deconv, modelGeneVar)


#add the batch column
#all.sce.deconv[[1]]$batch<-rep("ctrl_1",ncol(all.sce.deconv[[1]]))
all.sce.deconv[[1]]$batch<-rep("ctrl_2",ncol(all.sce.deconv[[1]]))
all.sce.deconv[[2]]$batch<-rep("isch",ncol(all.sce.deconv[[2]]))

  
colData(all.sce.deconv[[1]]) #check batch + label column

#take all the shared genes
allNames <- unlist(lapply(all.sce.deconv, function(x){rownames(x)}))
allNamesNb <- table(allNames)
universe <- names(allNamesNb)[allNamesNb==2]


#or use universe <- Reduce(intersect, lapply(all.sce.deconv, rownames))

# Subsetting the SingleCellExperiment object.
uni.sce <- lapply(all.sce.deconv, function(x){x[universe,]})


# Also subsetting the variance modelling results, for convenience.
uni.dec <- lapply(all.dec, function(x){x[universe,]})

rescaled <- multiBatchNorm(uni.sce, batch = NULL) #batch is ignored if multiple objects are present


#Test if we have NA feature
which(is.na(rownames(uni.dec[[2]])))

# compute average variance components across samples
combined.dec <- combineVar(uni.dec)

chosen.hvgs <- combined.dec$bio > 0
combined.dec$chosenHvg <- chosen.hvgs
table(chosen.hvgs) # 8507 genes for ctrl1 + ctrl2 /8045 genes for ctrl1 + ctrl2 + isch /8495 for ctrl_2 + isch 

metadata(uni.dec[[1]])

fit.pbmc_one <- metadata(uni.dec[[1]])
fit.pbmc_two <- metadata(uni.dec[[2]])

plot(fit.pbmc_two$mean, fit.pbmc_two$var, xlab="Moyenne de l'expression log",
     ylab="Variance de l'expression log")
curve(fit.pbmc_two$trend(x), col="dodgerblue", add=TRUE, lwd=2)

gene_var_ggplot<-uni.dec[[1]]
gene_var_gg<-uni.dec[[2]]

gene_var_ggplot %>% 
  # convert to tibble for ggplot
  as_tibble() %>% 
  # make the plot
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_line(aes(y = tech), colour = "dodgerblue", size = 1) +
  labs(x = "Moyenne de l'expression log", y = "Variance de l'expression log ")+

gene_var_gg %>% 
  # convert to tibble for ggplot
  as_tibble() %>% 
  # make the plot
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_line(aes(y = tech), colour = "dodgerblue", size = 1) +
  labs(x = "Moyenne de l'expression log", y = "Variance de l'expression log ")

par(mfrow=c(1,2))

#one<-rescaled[[1]]
two<-rescaled[[1]]
three<-rescaled[[2]]

four<-cbind(logcounts(rescaled[[1]]),logcounts(rescaled[[2]]))


col_four<-rbind(colData(two),colData(three))

row_four<-rowData(two)[1:4] #it should be the same ! we don't carry with us the mean , detected , gene_sparsity

uncorrected <- SingleCellExperiment(list(logcounts=four),colData=col_four,rowData=row_four)

#check 
all(rownames(combined.dec) == rownames(uncorrected))

rowData(uncorrected) <- cbind(rowData(uncorrected), combined.dec)

set.seed(0010101010)
uncorrected <- runPCA(uncorrected,
                      subset_row=chosen.hvgs,
                      BSPARAM=BiocSingular::RandomParam())

# build shared nearest-neighbour graph
snn.gr <- buildSNNGraph(uncorrected,use.dimred="PCA") #because we have many cells (after mixing batches) we tend to have a higher k !

# identify cluster with the walk trap method
clusters_louv <- igraph::cluster_louvain(snn.gr)$membership
clusters_walk<-igraph::cluster_walktrap(snn.gr)$membership
clusters_leid<-igraph::cluster_leiden(snn.gr)$membership

#The idea here is to see if the clusterization with uncorrected values will produce clusters that are either from isch or ctrl ,in that case
#it means that we have a batch effect which is not du to biological differences (biaised techniques)

table(clusters_louv)
table(clusters_walk)
table(clusters_leid)


#we these following plots we want to see if we have only-ctrl cluster or vice versa (batch effect diagnosis) 
# get number of cells for each {cluster, batch} pair
clusterTab <- data.frame("clusters" = clusters_leid, "batch" = uncorrected$batch)

uncorrected


percent.var <- attr(reducedDim(uncorrected), "percentVar")


percent.var=as.data.frame(cbind(percent.var,seq(1:50)))

ggplot(percent.var, aes(x=V2, y=percent.var)) + 
  geom_point()+
  geom_vline(xintercept = 10.5,linetype =1,size=0.3,color="blue")+
  geom_hline(yintercept = 0.9,linetype =2,size=0.3,color="red")+
  xlab("Composantes principales")+ylab("% de variance expliquée")

plotReducedDim(uncorrected, dimred="PCA", colour_by="leiden")
plotReducedDim(uncorrected, dimred="PCA", ncomponents=3,
               colour_by="leiden")


factor(clusters)

ClusterInfo <- clusterTab %>% 
  group_by(clusters, batch) %>%
  summarise(cells = n()) 

p1 <- ggplot(data=ClusterInfo, aes(x=clusters,y=cells, fill=batch)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  ggtitle("uncorrected, cell numbers") +
  theme(legend.text = element_text(size = 7))
p2 <- ggplot(data=clusterTab, aes(x=clusters, fill = batch)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  ggtitle("uncorrected, proportions") +
  theme(legend.text = element_text(size = 7))

plot_grid(p1, p2, ncol=1)

#PLOT TSNE
set.seed(1111001)
uncorrected <- runTSNE(uncorrected, dimred = "PCA")
uncorrected<-runUMAP(uncorrected,dimred="PCA")



# draw:
p.tsne <- plotTSNE(uncorrected_clone,
                   colour_by = "batch") +
  theme(legend.text = element_text(size = 7))

p.tsne

#Another
p.tsne + facet_wrap(~uncorrected$batch)


#Another brief one

plotUMAP(uncorrected,
         colour_by = "batch",text_by="label")


#get Some correction on the values
set.seed(1000101001)
mnn.out <- fastMNN(rescaled,
                   auto.merge = TRUE,
                   d = 50,
                   k = 20,
                   subset.row = chosen.hvgs , # to select all the genes for the analysis rep(TRUE,length(rownames(rescaled[[2]])))
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))


mnn.out$batch[which(mnn.out$batch==1)]<-rep("ctrl_2",length(which(mnn.out$batch==1)))
mnn.out$batch[which(mnn.out$batch==2)]<-rep("isch_1",length(which(mnn.out$batch==2)))

mnn.out$batch<-factor(mnn.out$batch)


mnn.out.corre.dim <- dim(reducedDim(mnn.out, "corrected"))
mnn.out.recon.dim <- dim(assay(mnn.out, "reconstructed"))

colDataList <- lapply(rescaled, function(x){colData(x)})
colDataDf <- do.call(rbind, colDataList)
colData(mnn.out) <- DataFrame(colDataDf)

#Mixing between batches



#we attribute the clusterization of mnn.out (after correction of the batch effects)to the uncorrected one for the marker gene analysis later 
snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected", k=20)
clust_walk <- igraph::cluster_walktrap(snn.gr)$membership
clust_louv<-igraph::cluster_louvain(snn.gr)$membership
clust_leid<-igraph::cluster_leiden(snn.gr)$membership





#louvain clust hesitate a lot it's dangerous (it gives 3 conformation (9c / 10c / 11c)

mnn.out$louvain<-factor(clust_louv) #with that linecode
mnn.out$walktrap<-factor(clust_walk)
mnn.out$leiden<-factor(clust_leid)

table(clust_louv)

table(clust_walk)

table(clust_leid)

?igraph::cluster_louvain
#Now it's time to replot the porportions chart of each clusters and see wether the ish-only clusters or ctrl-only clusters disapeared


clusterTab.mnn <- data.frame(clusters=clust_louv, batch=mnn.out$batch)

ClusterInfo.mnn <- clusterTab.mnn %>% 
  as.data.frame() %>%
  group_by(clusters,batch) %>%
  summarise(cells = n())

mp1 <- ggplot(data=ClusterInfo.mnn, aes(x=clusters,y=cells, fill=batch)) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  geom_col() +
  theme(legend.text = element_text(size = 7))
mp2 <- ggplot(data=clusterTab.mnn, aes(x=clusters, fill=batch)) +
  geom_bar(position = "fill") +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.text = element_text(size = 7))

plot_grid(mp1, mp2, ncol=1)

#how much variance did I lost here .
round(metadata(mnn.out)$merge.info$lost.var,2)

#we add the louvain partitionning after the MNN correction (it's the more precise)
colData(uncorrected)$louvain <- factor(colData(mnn.out)$louvain)

#we pass the walktrap clustering to the uncorrected one ( when i got the isch patient i decided to not carry the walktrap anymore)
colData(uncorrected)$walktrap<-factor(clust_walk)

#same for leiden method
colData(uncorrected)$leiden<-factor(clust_leid)

# Marker gene identification (biased analysis with pvals , it's better to follow OSCA workflow ) ----

# identify marker genes based on mean expression differences
# default options do not need to be specified, but shown here for illustration
markers_default <- findMarkers(
  uncorrected, 
  groups = factor(uncorrected$louvain), # clusters to compare
  block = uncorrected$batch,    # covariates in statistical model
  test.type = "t",   # t-test (default)
  direction = "any", # test for either higher or lower expression (default)
  lfc = 0, # null hypothesis log-fold-change = 0 (default)
  pval.type = "any" # ranking of p-values based on any comparison (default)
) 

# returns a list of length equal to the number of clusters
markers_default

# check the result of a particular cluster
rownames(markers_default[[10]][markers_default[[10]]$Top<=1,])

markers_default[[1]]


# Others ----


set.seed(155478811) #in order to have the same conformation in TSNE
mnn.out <- runTSNE(mnn.out, dimred="corrected")
mnn.out<-runUMAP(mnn.out,dimred="corrected")

par(mfrow=c(1,2))


table(uncorrected$batch)


plotUMAP(mnn.out, 
         colour_by = "batch",
         text_by = "louvain",
         by_exprs_values = "reconstructed")

plotTSNE(mnn.out, 
         colour_by = "louvain",
         by_exprs_values = "reconstructed")


table(clusters.man)

gridExtra::grid.arrange(
  plotTSNE(mnn.out, colour_by = "Tnnt2",by_exprs_values = "reconstructed"),
  plotTSNE(mnn.out, colour_by = "louvain",by_exprs_values = "reconstructed"),
  ncol=2
)

table(mnn.out$label)

?grid.arrange
?plotTSNE
table(mnn.out$label)

which(rownames(mnn.out)=="MRF1")

markers_default[[10]][which(rownames(markers_default[[10]])=="PAX7"),]

which(rownames(markers_default[[10]]["GNLY",]))

which(rownames(markers_default[[10]])=="PAX7")



##########Research of DE genes between conditions###############


#i want to count matrix to uncorrected in order to follow the OSCA DE analysis
assay(uncorrected,"counts")<-cbind(counts(rescaled[[1]]),counts(rescaled[[2]]))


#subset a sce which contains only satellite cells
sce.subset_endo<-uncorrected[,which(uncorrected$leiden=="4")]
sce.subset_pcv_endo<-uncorrected[,which(uncorrected$leiden=="10")]
sce.subset_sat<-uncorrected[,which(uncorrected$leiden=="7")]


#we are going to work as though we have pseudo-bulk samples!!!
summed<-aggregateAcrossCells(sce.subset, 
                             id=colData(sce.subset)[,c("leiden","batch")])



library(edgeR)

igrec <- DGEList(counts(summed), samples=colData(summed)[,c(2,3,4,10,12)])
igrec

?DGEList

keep <- filterByExpr(igrec, group=igrec$samples$batch)

#i can try to work with all the genes , because i can aford it computationally
#i can also try to work with subset provided by the variable keep
#igrec<-igrec[keep,]

igrec <- calcNormFactors(igrec)
igrec

par(mfrow=c(1,2))
for (i in seq_len(ncol(igrec))) {
  plotMD(igrec, column=i)
}

plotMDS(cpm(igrec, log=TRUE),col=ifelse(igrec$samples$batch, "red", "blue"))
#we can use MDS if we have more than 2D.

#just for visualization
plot(cpm(igrec,log = TRUE))


rownames(igrec$samples)<-c("clust_4_ctrl","clust_4_isch","clust_7_ctrl","clust_7_isch","clust_10_ctrl","clust_10_isch")
colnames(igrec)<-c("clust_4_ctrl","clust_4_isch","clust_7_ctrl","clust_7_isch","clust_10_ctrl","clust_10_isch")

design <- model.matrix(~factor(batch), igrec$samples)
design


#Estimate the NB dispersion ( here we want to model mean-variance trend) 
igrec <- estimateDisp(igrec, design)
summary(igrec$trended.dispersion)

plotBCV(igrec)

#estimate quasi-likelihood
fit <- glmQLFit(igrec, design, robust=TRUE, abundance.trend=FALSE)
summary(fit$var.prior)

plotQLDisp(fit,cex=1)

?glmQLFTest

res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))
#no DE foud so  i have a problem


####use directly pseudobulkDGE#######

igrec$samples$batch<-c("endo","sat")
de.results <- pseudoBulkDGE(summed, 
                            label=c("4","7"),
                            design=~factor(batch)
                            #coef="ctrl_2",
                            #condition=summed.filt$tomato 
)

#no result here so i have to test other things


###Compare clusters with scoreMarkers and/or FindMarkers###########



trip.sce<-all.sce.deconv[[2]]



marker_news<-scoreMarkers(sce.subset_endo,sce.subset_endo$batch)
marker_news<-scoreMarkers(trip.sce,trip.sce$leiden)


marker_news


khayer <- marker_news[[9]]
ordered <- khayer[order(khayer$mean.AUC, decreasing=TRUE),]
sdd<-as.data.frame(cbind(ordered,rownames(ordered))) # showing basic stats only, for brevity.
#ordered[1:10,]
View(sdd)

table(trip.sce$leiden)



sce.subset_endo<-uncorrected[,which(uncorrected$leiden=="4")]

sce.subset_faps<-uncorrected[,which(uncorrected$leiden=="3")]

markers_default <- findMarkers(
  sce.subset_faps, 
  groups = factor(sce.subset_faps$batch), # clusters to compare
  test.type = "t",   # t-test (default)
  direction = "any", # test for either higher or lower expression (default)
  lfc = 0, # null hypothesis log-fold-change = 0 (default)
  pval.type = "any" # ranking of p-values based on any comparison (default)
)



# returns a list of length equal to the number of clusters
markers_default

View(cbind(rownames(markers_default[[2]]),as.data.frame(markers_default[[2]])))

#the first [[1]] is ctrl-isch and the second [[2]] is isch-ctrl (the right one for us)

inspect_sat<-cbind(as.data.frame(markers_default[[1]]),rownames(markers_default[[1]]))
#i used the cluster louvain to get the sat cells (leiden finds 222 and louvain finds 223)
#i don't select the markers_default[[2]] because it's the same as the one

colnames(inspect_sat)[6]<-"gene_name"


inspect_endo<-cbind(as.data.frame(markers_default[[2]]),rownames(markers_default[[2]]))
colnames(inspect_endo)[6]<-"gene_name"





