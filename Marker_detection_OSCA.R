##OSCA Book 


marker.info <- scoreMarkers(uncorrected,uncorrected$louvain)

marker.info

colnames(marker.info[["3"]]) # these are all the effect size summaries



chosen <- marker.info[["5"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
sdd<-as.data.frame(cbind(ordered[,1:4],rownames(ordered))) # showing basic stats only, for brevity.
ordered[1:10,1:4]

which(rownames(ordered)=="FLT4")

plotExpression(uncorrected, features=head(rownames(ordered)), 
               x="louvain", colour_by="louvain")

which(rownames(cohen.only)=="ACTA2")



write.csv(as.data.frame(ordered[0:100,1:4]),"cluster_9.csv",row.names = TRUE)
          
##Cohen's Distance

cohen.only <- chosen[,grepl("logFC.cohen", colnames(chosen))]
cohen.only<-cohen.only[order(cohen.only$mean.logFC.cohen,decreasing=TRUE),]
cohen.only[1:11,]

##LogFC

detect.only <- chosen[,grepl("logFC.detected", colnames(chosen))]
detect.only<-detect.only[order(detect.only$mean.logFC.detected,decreasing=TRUE),]
detect.only[1:11,]



#Save the first 20 markers 

setwd("/home/izem/Desktop/STAGE/Markers_Clusters/")

p= 20 #how many markers do you want to save

for(i in 1:length(marker.info)){
  chosen <- marker.info[[i]]
  ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
  stock<-as.data.frame(cbind(rownames(ordered)[1:p],ordered[1:p,]))
  write_csv(stock,paste("AUC_cluster_",i,".csv",sep = "")) # showing basic stats only, for brevity.
  cohen.only <- chosen[,grepl("logFC.cohen", colnames(chosen))]
  cohen.only<-cohen.only[order(cohen.only$mean.logFC.cohen,decreasing=TRUE),]
  stock<-as.data.frame(cbind(rownames(cohen.only)[1:10],cohen.only[1:p,]))
  write_csv(stock,paste("Cohen_cluster_",i,".csv",sep = ""))
  detect.only <- chosen[,grepl("logFC.detected", colnames(chosen))]
  detect.only<-detect.only[order(detect.only$mean.logFC.detected,decreasing=TRUE),]
  stock<-as.data.frame(cbind(rownames(ordered)[1:10],detect.only[1:p,]))
  write_csv(stock,paste("logFC_cluster_",i,".csv",sep = ""))
}


write_csv(as.data.frame(ordered[1:10]),paste("AUC_cluster_",i,".csv",sep = ""),row.names = TRUE)


cbind(rownames(ordered)[1:10],ordered[1:10,])



idea<-as.data.frame(cbind(uncorrected$batch,uncorrected$louvain))


#####Do this for SCSA (python)########

#scran_pbmc_3k.csv was generated by following command from sce object(due to its pairwise comparisons, we use the mean LFC instead):


markers <- findMarkers(uncorrected, uncorrected$louvain, pval.type="all")
res <- data.frame()
for (i in names(markers)){
  predata <- subset(markers[[i]],select=c(p.value,FDR))
  meandata <- as.matrix(apply(subset(markers[[i]],select=-c(p.value,FDR)),1,mean)) 
  if (length(res) == 0){
    colnames(meandata) <- paste("LFC",i,sep="_")
    colnames(predata) <- paste(names(predata),i,sep="_")
    res <- cbind(predata,meandata)
  }else{
    predata <- predata[rownames(res),]
    meandata <- as.matrix(meandata[rownames(res),])
    colnames(meandata) <- paste("LFC",i,sep="_")
    colnames(predata) <- paste(names(predata),i,sep="_")
    res <- cbind(res,predata,meandata)
  }
}

write.csv(res,file="/home/izem/Desktop/STAGE/SCSA/findMarkers.csv",quote=FALSE)



markers[[4]]
marker.info[[4]]



#Try scCATCH############

#il catch rien le fr??re

library(scCATCH)

WWE<-createscCATCH(logcounts(uncorrected),as.character(uncorrected$louvain))

findmarkergene(
  WWE,
  species = "Human",
  cluster = "9",
  if_use_custom_marker = FALSE,
  marker = cellmatch,
  cancer = "Normal",
  tissue = "Skeletal muscle",
  use_method = "1",
  comp_cluster = 6,
  cell_min_pct = 0.25,
  logfc = 0.25,
  pvalue = 0.05,
  verbose = TRUE
)

#####Try SingleR#######
library(celldex)
ref.data <- HumanPrimaryCellAtlasData(ensembl=FALSE)



library(SingleR)
predictions <- SingleR(test=uncorrected, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main,clusters = uncorrected$louvain)

chafou<-as.data.frame(predictions)

##From Genes to GO term########

#first get the vector of the 20 top genes from cluster 9

AUC_nine<-rownames(ordered)[1:20]

#Cohen_nine<-rownames(cohen.only)[1:20]
#logfc_nine<-rownames(detect.only)

tabla<-as.data.frame(toTable(org.Hs.egSYMBOL2EG[AUC_nine]))

maida<-as.data.frame(toTable(org.Hs.egGO))


library(org.Hs.eg.db)

GOGOGO<-merge(tabla, maida, by.x=
        "gene_id", by.y="gene_id")


#Stock into 3 dataframes.
nine_MF<-GOGOGO[which(GOGOGO$Ontology=="MF"),]

nine_CC<-GOGOGO[which(GOGOGO$Ontology=="CC"),]

nine_BP<-GOGOGO[which(GOGOGO$Ontology=="BP"),]


#find the explanation of go term 

library(GO.db)


#Term("GO:0016310")
#Definition("GO:0016310")

?Definition
gros_bucket<-list(nine_BP,nine_CC,nine_MF)

for(n in 1:3){
  current=NULL
  terme=NULL
  definition=NULL
  current<-gros_bucket[[n]]
  go_idea<-current$go_id
  for (m in 1:length(current$go_id)){
    terme[m]<-Term(go_idea[m])
    definition[m]<-Definition(go_idea[m])
  }
  current<-cbind(current,terme,definition)
  gros_bucket[[n]]<-current
  }


nine_BP<-gros_bucket[[1]]
nine_CC<-gros_bucket[[2]]
nine_MF<-gros_bucket[[3]]

setwd("/home/izem/Desktop/STAGE/Cluster_9/")
write_csv(nine_BP,paste("BP_9",".csv",sep = ""))
write_csv(nine_CC,paste("CC_9",".csv",sep = ""))
write_csv(nine_MF,paste("MF_9",".csv",sep = ""))





###Investigations ##########

#let's check how many Faps cells do we have in isch 

sub_idea<-idea[which(idea$V2==1),]

length(sub_idea$V1)

propor_ctrl_1<-NULL
propor_ctrl_2<-NULL

for(i in 1:14){
  sub_idea<-idea[which(idea$V2==i),]
  propor_ctrl_1[i]<-(length(sub_idea[which(sub_idea$V1=="isch"),]$V1)/length(sub_idea$V1))*100
  propor_ctrl_2[i]<-(length(sub_idea[which(sub_idea$V1=="ctrl"),]$V1)/length(sub_idea$V1))*100
}

final_propor<-data.frame("percentage_ctrl_1"=propor_ctrl_1,"percentage_ctrl_2"=propor_ctrl_2)


