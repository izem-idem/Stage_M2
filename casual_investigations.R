

###let's wich of the clustering is the best


table(clust_louv)

table(clust_walk)

table(clust_leid)


plotUMAP(mnn.out, 
         colour_by = "PAX7",
         text_by = "louvain",
         by_exprs_values = "reconstructed")

gridExtra::grid.arrange(
  plotTSNE(mnn.out, colour_by = "CDKN1A",text_by = "louvain",by_exprs_values = "reconstructed"),
  plotTSNE(mnn.out, colour_by = "batch",text_by = "leiden",by_exprs_values = "reconstructed"),
  plotTSNE(mnn.out, colour_by = "batch",text_by = "walktrap",by_exprs_values = "reconstructed"),
  ncol=2
)


gridExtra::grid.arrange(
  plotTSNE(test, colour_by = "leiden",by_exprs_values = "reconstructed"),
  plotTSNE(test, colour_by = "batch",text_by = "leiden",by_exprs_values = "reconstructed"),
  ncol=2
)

test<-mnn.out


test$batch<-as.factor(test$batch)

levels(test$batch)<-c("CTRL","IMI")


?aggregate


###Check markers#####

marker.info <- scoreMarkers(uncorrected,uncorrected$leiden)

marker.info

colnames(marker.info[["3"]]) # these are all the effect size summaries

chosen <- marker.info[["4"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
#sdd<-as.data.frame(cbind(ordered[,1:4],rownames(ordered))) # showing basic stats only, for brevity.
ordered[1:10,1:4]

which(rownames(ordered)=="FLT1")




library(celldex)
ref.data <- HumanPrimaryCellAtlasData(ensembl=FALSE)


library(SingleR)
predictions <- SingleR(test=uncorrected, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main,clusters = uncorrected$leiden)

chafou_leiden<-as.data.frame(predictions)




###Go-term search#########

library(org.Hs.eg.db)
library(GO.db)


for(i in c(4,7,10)){
  chosen <- marker.info[[i]]
  ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
  douka<-rownames(ordered)[1:20] #get the 20th first markers
  tabla<-as.data.frame(toTable(org.Hs.egSYMBOL2EG[douka]))
  maida<-as.data.frame(toTable(org.Hs.egGO))
  GOGOGO<-merge(tabla, maida, by.x=
                  "gene_id", by.y="gene_id")
  #stock it into 3 different dataframes
  nine_MF<-GOGOGO[which(GOGOGO$Ontology=="MF"),]
  nine_CC<-GOGOGO[which(GOGOGO$Ontology=="CC"),]
  nine_BP<-GOGOGO[which(GOGOGO$Ontology=="BP"),]
  
  
}


#first get the vector of the 20 top genes from cluster 9

AUC_nine<-rownames(ordered)[1:20]

#Cohen_nine<-rownames(cohen.only)[1:20]
#logfc_nine<-rownames(detect.only)

tabla<-as.data.frame(toTable(org.Hs.egSYMBOL2EG[AUC_nine]))

maida<-as.data.frame(toTable(org.Hs.egGO))




GOGOGO<-merge(tabla, maida, by.x=
                "gene_id", by.y="gene_id")


#Stock into 3 dataframes.
nine_MF<-GOGOGO[which(GOGOGO$Ontology=="MF"),]

nine_CC<-GOGOGO[which(GOGOGO$Ontology=="CC"),]

nine_BP<-GOGOGO[which(GOGOGO$Ontology=="BP"),]


#find the explanation of go term 



#Term("GO:0016310")
#Definition("GO:0016310")

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


#####Rearrange GO_DF ##########

col_up=list()
col_down=list()

library(stringr) #to use str_trim

for(i in 1:nrow(endo_GO_df)){ 
  vect_col_up=NULL
  vect_col_down=NULL
  x=0
  ester=str_split(endo_GO_df[i,8],",") #get the gene_id column and split it into a vector
  for(j in 1:length(ester[[1]])){
    x=x+1
    gene=str_trim(ester[[1]][j],"both") #remove spaces at start and end of the string
    if(circ_input_two[circ_input_two$ID==gene,2]>0){ #if gene is up regulated then 
      vect_col_up[j]=gene
    }
    if(circ_input_two[circ_input_two$ID==gene,2]<0){
      vect_col_down[j]=gene
    }
  }
  col_up[[i]]=vect_col_up
  col_down[[i]]=vect_col_down
}


GO_clone=GO_df[,c(1,2,3,6)]

GO_clone$up=rep(NA,nrow(GO_clone))
GO_clone$down=rep(NA,nrow(GO_clone))

for(i in 1:nrow(GO_clone)){
  GO_clone[i,5]=toString(na.omit(col_up[[i]]))
  GO_clone[i,6]=toString(na.omit(col_down[[i]]))
}

write.csv(GO_clone,"Cbon.csv")


#### Working just with entrez ID because of satané mitochondriale genes


col_up=list()
col_down=list()

library(stringr) #to use str_trim

for(i in 1:nrow(endo_GO_df)){ 
  vect_col_up=NULL
  vect_col_down=NULL
  x=0
  ester=str_split(endo_GO_df[i,8],"/") #get the gene_id column and split it into a vector
  for(j in 1:length(ester[[1]])){
    x=x+1
    gene=str_trim(ester[[1]][j],"both") #remove spaces at start and end of the string
    if(lion[lion$ENTREZID==gene,7]>0){ #if gene is up regulated then 
      vect_col_down[j]=gene
    }
    if(lion[lion$ENTREZID==gene,7]<0){
      vect_col_up[j]=gene
    }
  }
  col_up[[i]]=vect_col_up
  col_down[[i]]=vect_col_down
}


GO_clone=endo_GO_df[,c(1,2,3,6)]

GO_clone$up=rep(NA,nrow(GO_clone))
GO_clone$down=rep(NA,nrow(GO_clone))

for(i in 1:nrow(GO_clone)){
  GO_clone[i,5]=toString(na.omit(col_up[[i]]))
  GO_clone[i,6]=toString(na.omit(col_down[[i]]))
}

for(i in 1:nrow(GO_clone)){
  GO_clone[i,5]=toString(bitr(strsplit(gsub(" ","",GO_clone[i,5]),",")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL)
  GO_clone[i,6]=toString(bitr(strsplit(gsub(" ","",GO_clone[i,6]),",")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL)
}

toString(bitr(strsplit(gsub(" ","",GO_clone[1,6]),",")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL)

