rm(list = ls())

library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)

setwd("/home/izem/Desktop/STAGE/")


bp.params <- MulticoreParam(workers = 3)

sample.path <- "/home/izem/Desktop/STAGE/Izem_Single-cell/muscle_1_isch/outs/raw_feature_bc_matrix/"

sample.path.two <- "/home/izem/Desktop/Master/STAGE_M2/raw_feature_bc_matrix_muscle_2_ctrl/"

sample.path.test <- "/home/izem/Desktop/Master/STAGE_M2/count-mkref_Explore/run_count_1kpbmcs(with_index)/outs/raw_feature_bc_matrix/"
  
sce.sing <- read10xCounts(sample.path, col.names=TRUE, BPPARAM = bp.params)

sce.sing.two <- read10xCounts(sample.path.two, col.names=TRUE, BPPARAM = bp.params)

sce.sing #there is a many genes that are not expressed in any cells , we can remove those.

sce.sing.two # we can see that cellranger count puts it at 778 cells


rowData(sce.sing) #to acces information about features

colData(sce.sing) #to access information about barcodes.


####Number of genes detected per cell

genesPerCell <- colSums(counts(sce.sing) > 0)

genesPerCell.two <- colSums(counts(sce.sing.two) > 0)

genesPerCell
summary(genesPerCell)
#on a une médiane de 1 géne exprimé pour les 414 106 eventuelles barcodes
#bon devant autant de barcodes je comprends que le signal est dilué
#on va voir comment ça évolue. 

summary(genesPerCell.two)

plot(density(genesPerCell), main="", xlab="Genes per cell")

######Total UMI for a gene versus the number of times detected  
#on voit que des génes qui sont fortement exprimés ont tendance à être présent dans les cellules.


#Counts stocke les informations dans une matrice sparse ie qui contient beaucoup de zéros 
#le % de zéro dans cette matrice c'est le % de sparsité et le % de non-zéro c'est le % de densité

tmpCounts <- counts(sce.sing)[,1:1000]

plot(rowSums(tmpCounts),
     rowMeans(tmpCounts > 0),
     log = "x",
     xlab="total number of UMIs",
     ylab="proportion of cells expressing the gene"
)
rm(tmpCounts)


###### On va essayer d'enlever quelques features pour alléger la matrice sinon ça va être problématiques (SAHMRI Workshop)

set.seed(1556)
e.out <- emptyDrops(counts(sce.sing),BPPARAM=bp.params)
sce.pbmc <- sce.sing[,which(e.out$FDR <= 0.005)] #on récupére les indices des barcodes qui ont un FDR inférieur ou égal à 0.001
unfiltered <- sce.pbmc
sce.pbmc

#donc là j'élimine toutes les cellules qui sont considérés comme vide.

####Remove Undetected genes (BROAD)

detected_genes <- rowSums(counts(sce.pbmc)) > 0
table(detected_genes)

sce.pbmc <- sce.pbmc[detected_genes,]



######Distribution of counts for a gene across cells (BROAD) ####


rel_expression <- t( t(counts(sce.pbmc)) / colSums(counts(sce.pbmc))) * 100
rownames(rel_expression) <- rowData(sce.pbmc)$Symbol
most_expressed <- sort(rowSums( rel_expression ),T)
plot_data <- as.matrix(t(rel_expression[names(most_expressed),]))

write.table(head(t(plot_data),n=10L), file="file.csv", row.names = T, col.names = T)
#je fais la transposé car je veux avoir en col les cells et en row les genes (sinon le fichier excel ne s'affiche pas en entier)
#ici je choisis d'afficher les 10e lignes de la matrice dans le fichier excel (mais je peux faire plus bien sûr)


boxplot(head(plot_data,n=10L), cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE)
#je n'arrive pas à faire ce plot du fait de ça grande taille




#####Modify the droplet annotation

#we'll comeback to this part when we'll have many samples

#######Anottate genes (BROAD) ####### 

ah <- AnnotationHub()
ens.mm.98 <- query(ah, c("Homo sapiens", "EnsDb", 98))[[1]] 

genes <- rowData(sce.pbmc)$ID
gene_annot <- AnnotationDbi::select(ens.mm.98, 
                                    keys = genes,
                                    keytype = "GENEID",
                                    columns = c("GENEID", "SEQNAME"))
rowData(sce.pbmc) <- merge(rowData(sce.pbmc), gene_annot, by = "ID", sort=FALSE)
rownames(rowData(sce.pbmc)) <- rowData(sce.pbmc)$ID




#####Add per cell QC metrics (BROAD) ######

is.mito <- which(rowData(sce.pbmc)$Chromosome=="MT")
sce.pbmc<- addPerCellQC(sce.pbmc, subsets=list(Mito=is.mito), BPPARAM = bp.params)


colData(sce.pbmc)
########QC_metrics distribution########


###For multiple samples##
plotColData(sce, x="Sample", y="sum",other_fields="SampleGroup") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Total count")

plotColData(sce, x="Sample", y="detected", other_fields="SampleGroup") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Detected features")

plotColData(sce, x="Sample", y="subsets_Mito_percent", other_fields="SampleGroup") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") +
  ggtitle("Mito percent")


colData(sce.pbmc) %>% 
  as.data.frame() %>% 
  arrange(subsets_Mito_percent) %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(colour = subsets_Mito_percent > 10)) + 
  facet_wrap(vars(SampleGroup))


####For one sample###


plotColData(sce.pbmc, y="sum") +
  scale_y_log10() + ggtitle("Total count")

plotColData(sce.pbmc, x="sum", y="subsets_Mito_percent")

gridExtra::grid.arrange(
  plotColData(sce.pbmc, y = "sum") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce.pbmc, y = "detected") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce.pbmc, y = "subsets_Mito_percent") + ggtitle("Mito percent"),
  ncol = 2
)


#########Identification of low-quality cells (BROAD)########

###############################################With library size ###

qqnorm(colData(sce.pbmc)$sum)
qqline(colData(sce.pbmc)$sum,col=2)
#FIRST CHECK THE NORMALITY OF THIS DATA , BECAUSE THE FOLLOWING Method (isoutlier)
#assumes that we have a gausssian distribution

shapiro.test(colData(sce.pbmc)$detected)

#subsets_Mito_percent
#Even if the data are not normally distributed , if the distribution is 
#symmetric we can apply the MAD based approach
library(moments)

skewness(colData(sce.pbmc)$subsets_Mito_percent) 

#this score should be the closest to 0 .

low_lib_size <- isOutlier(sce.pbmc$sum, log=TRUE, type="lower")
table(low_lib_size)

?isOutlier
attr(low_lib_size, "thresholds")
#we can see the threshold applied to define a cell as low quality one

colData(sce.pbmc)$low_lib_size <- low_lib_size

##############For one sample
plotColData(sce.pbmc, 
            y="sum",
            colour_by = "low_lib_size") +
  scale_y_log10() + 
  labs(y = "Total count", title = "Total count") +
  guides(colour=guide_legend(title="Discarded"))

##############For multiple sample
colData(sce.pbmc)$low_lib_size <- low_lib_size
plotColData(sce.pbmc, 
            x="Sample", 
            y="sum",
            other_fields="SampleGroup", 
            colour_by = "low_lib_size") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  labs(y = "Total count", title = "Total count") +
  guides(colour=guide_legend(title="Discarded"))


############################################### Number of genes ### 

low_n_features <- isOutlier(sce.pbmc$detected, log=TRUE, type="lower")
table(low_n_features)

attr(low_n_features, "thresholds")[1]

#we remove low quality cells that do not express more than x genes

colData(sce.pbmc)$low_n_features <- low_n_features


##############For one sample
plotColData(sce.pbmc, 
            y="detected",
            colour_by = "low_n_features") + 
  scale_y_log10() + 
  labs(y = "Genes detected", title = "Genes detected") +
  guides(colour=guide_legend(title="Discarded"))


##############For multiple sample
plotColData(sce.pbmc, 
            x="Sample", 
            y="detected",
            other_fields="SampleGroup", 
            colour_by = "low_n_features") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  labs(y = "Genes detected", title = "Genes detected") +
  guides(colour=guide_legend(title="Discarded"))

###############################################Mitochondrial content

high_Mito_percent <- isOutlier(sce.pbmc$subsets_Mito_percent, type="higher")
table(high_Mito_percent)

attr(high_Mito_percent, "thresholds")[2]

colData(sce.pbmc)$high_Mito_percent <- high_Mito_percent

##############For one sample

plotColData(sce.pbmc,  
            y="subsets_Mito_percent",
            colour_by = "high_Mito_percent") + 
  labs(y = "Percentage mitochondrial UMIs",
       title = "Mitochondrial UMIs") +
  guides(colour=guide_legend(title="Discarded"))



##############For multiple sample

plotColData(sce.pbmc,  
            x="Sample",
            y="subsets_Mito_percent",
            other_fields="SampleGroup",
            colour_by = "high_Mito_percent") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  labs(y = "Percentage mitochondrial UMIs",
       title = "Mitochondrial UMIs") +
  guides(colour=guide_legend(title="Discarded"))


#To see how much low quality cells did we removed
data.frame(`Library Size` = sum(low_lib_size),
           `Genes detected` = sum(low_n_features),
           `Mitochondrial UMIs` = sum(high_Mito_percent),
           Total = sum(low_lib_size | low_n_features | high_Mito_percent))


sce.pbmc

###################TO DO IT ALL AT ONCE
#we put subset_mito_percent to add it to the QC results
cell_qc_results <- quickPerCellQC(colData(sce.pbmc),percent_subsets=c("subsets_Mito_percent"))

colSums(as.data.frame(cell_qc_results))

sce.pbmc$discard <- cell_qc_results$discard



#################To remove low Q cells 

sce.filtered <- sce.pbmc[, !sce.pbmc$discard]

#We need to remove ancient (pre_filtered) QC metrics and recalculate them

colData(sce.filtered) <- colData(sce.filtered)[,1:3] # remove all QC col and keep the first 3 cols

sce.filtered <- addPerCellQC(sce.filtered, BPPARAM = bp.params)







##############QC_normalized##########

colData(sce.pbmc)$log10_counts<-log10(sce.pbmc$sum +1)

low.lib <- isOutlier(sce.pbmc$log10_counts, nmads=3, type="lower")

colData(sce.pbmc)$log10_detected<-log10(sce.pbmc$detected +1)

low.nfeat <- isOutlier(sce.pbmc$log10_detected, nmads=3, type="lower")


discard <- low.lib | low.nfeat
data.frame(LowLib=sum(low.lib), LowNFeat=sum(low.nfeat),discard=sum(discard))

######################QC_based on sparsity (BROAD)######


#we have to remember that "sum" and "detected" are log-transformed (when we want to  compare them)


#We'll search for : 
##the cell sparsity: for each cell, the proportion of genes that are not detected
## the gene sparsity: for each gene, the proportion of cells in which it is not detected

sce.pbmc <- addPerFeatureQC(sce.pbmc, BPPARAM = bp.params)
rowData(sce.pbmc)

# mean - the mean UMI count for the gene across all cells
# detected - the percentage of cells in which the gene was detected

head(sort(rowData(sce.pbmc)[,5],decreasing = T),n=600L) # best 100 genes which have the highest mean

head(sort(rowData(sce.pbmc)[,6],decreasing = T),n=600L) # most 100 detected genes

colData(sce.pbmc)

#we can use sce.pbmc or sce.filtered (we can try both)


colData(sce.pbmc)$cell_sparsity <- 1 - (colData(sce.pbmc)$detected / nrow(sce.pbmc))
rowData(sce.pbmc)$gene_sparsity <- (100 - rowData(sce.pbmc)$detected) / 100




hist(sce.pbmc$cell_sparsity, breaks=50, col="grey80", xlab="Cell sparsity", main="")

hist(rowData(sce.pbmc)$gene_sparsity, breaks=50, col="grey80", xlab="Gene sparsity", main="")


#####FILTER BY SPARSITY (BROAD)#####

sparse.cells <- sce.pbmc$cell_sparsity > 0.995
mito.cells <- sce.pbmc$subsets_Mito_percent >= 10

min.cells <- 1 - (20 / ncol(sce.pbmc))
sparse.genes <- rowData(sce.pbmc)$gene_sparsity > min.cells



table(sparse.genes)

table(sparse.cells, mito.cells)

table(mito.cells)
table(sparse.cells)

sce.pbmc$cell_sparsity



