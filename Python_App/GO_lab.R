library(org.Hs.eg.db)
library(GO.db)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(enrichplot)




args <- commandArgs(TRUE)
print(args[1])
setwd(args[1])

#load("/home/izem/Desktop/STAGE/bouge_pas_yaw.RData")


print(args[2])
if(args[2]=="GO"){ #meaning GO
  print("GO Module")
  filou<-read.delim2("target",header = FALSE)
  print(args[6])
  if(args[6]=="human"){orga="org.Hs.eg.db";gene<-read.csv("Human_gene.csv")}
  else{orga="org.Mm.eg.db";gene<-read.csv("Mouse_gene.csv")}

  filou<-merge(filou,gene,by.x="V1",by.y="Symbol")
  #print(one)
  two<-bitr(filou$ID,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = orga)
  print(args[4])
  if(args[4]=="1"){
    print("Background gene list detected for GO")
    fila<-read.delim2("background",header = FALSE) #here we have a lot of genes no need to call human_gene or mouse_gene
    fila<-as.vector(t(fila))
    print("flag three")
    three<-bitr(fila,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = orga)$ENSEMBL
    print("flag_four")
    four<-bitr(three,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = orga)
    print("flag_chebab")
    chebab_GO<-enrichGO(two$ENTREZID,OrgDb = orga,universe =four$ENTREZID,pvalueCutoff = as.numeric(args[7]),pAdjustMethod =args[8],ont =args[9])
    print(chebab_GO)
    chebab_GO_df<-as.data.frame(chebab_GO)
  }
  else{
    print("No Background gene list detected for GO")
    chebab_GO<-enrichGO(two$ENTREZID,OrgDb = orga,pvalueCutoff =as.numeric(args[7]),pAdjustMethod = args[8],ont = args[9])
    print(chebab_GO)
    chebab_GO_df<-as.data.frame(chebab_GO)
  }

  for(i in 1:nrow(chebab_GO_df)){
    chebab_GO_df[i,8]=toString((bitr(strsplit(chebab_GO_df[i,8],"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = orga)$SYMBOL))
    chebab_GO_df[i,2]=paste(chebab_GO_df[i,2]," (",chebab_GO_df[i,1],")")
  }


  write.csv(chebab_GO_df,"Enrich_result_GO.csv")

}else{ #meaning KEGG
  print("KEGG module")
  filou<-read.delim2("target",header = FALSE)

  if(args[6]=="human"){orga="org.Hs.eg.db";org="hsa";gene<-read.csv("Human_gene.csv")}#i use this one because the transition is more acurate
  #bitr function loses one or more genes in the translation
  else{orga="org.Mm.eg.db";org="mmu";gene<-read.csv("Mouse_gene.csv")}

  filou<-merge(filou,gene,by.x="V1",by.y="Symbol")
  two<-bitr(filou$ID,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = orga)


  if(args[4]==1){
    print("Background gene list detected for KEGG")
    fila<-read.delim2("background",header = FALSE)
    fila<-as.vector(t(fila))
    three<-bitr(fila,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = orga)$ENSEMBL
    four<-bitr(three,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = orga)
    chebab_KEGG<-enrichKEGG(two$ENTREZID,organism = org,universe =four$ENTREZID,pvalueCutoff = as.numeric(args[7]),pAdjustMethod =args[8],keyType="ncbi-geneid")
    print(chebab_KEGG)
    chebab_KEGG_df<-as.data.frame(chebab_KEGG)
  }
  else{
    print("No Background gene list detected for KEGG")
    chebab_KEGG<-enrichKEGG(two$ENTREZID,organism = org,pvalueCutoff = as.numeric(args[7]),pAdjustMethod = args[8],keyType="ncbi-geneid")
    print(chebab_KEGG)
    chebab_KEGG_df<-as.data.frame(chebab_KEGG)
  }

  for(i in 1:nrow(chebab_KEGG_df)){
    chebab_KEGG_df[i,8]=toString((bitr(strsplit(chebab_KEGG_df[i,8],"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = orga)$SYMBOL))
    chebab_KEGG_df[i,2]=paste(chebab_KEGG_df[i,2]," (",chebab_KEGG_df[i,1],")")
  }
  write.csv(chebab_KEGG_df,"Enrich_result_KEGG.csv")
}
