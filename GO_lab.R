library(org.Hs.eg.db)
library(GO.db)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)

Trim_Print <- function(endo_GO_df,prenom) {
  col_up=list()
  col_down=list()
  
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
  
  write.csv(GO_clone,paste(prenom,".csv"))
}	


  

args <- commandArgs(TRUE)


setwd(args[1])

if(args[2]==1){ #meaning GO
  
  filou<-read.delim2("text",header = FALSE)
  filou<-as.vector(t(filou))
  
  if(args[6]=="human"){orga="org.Hs.eg.db"}
  else{orga="org.Mm.eg.db"}

  one<-bitr(filou,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = orga)$ENSEMBL
  two<-bitr(one,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = orga )
  
  if(args[4]==1){
    fila<-read.delim2("background",header = FALSE)
    fila<-as.vector(t(fila))
    three<-bitr(fila,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = orga)$ENSEMBL
    four<-bitr(three,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = orga)
    chebab_GO<-enrichGO(two$ENTREZID,OrgDb = orga,universe =four$ENTREZID,pvalueCutoff = args[7],pAdjustMethod =args[8],ont =args[9])
    chebab_GO_df<-as.data.frame(chebab_GO)
  }

  else{
    chebab_GO<-enrichGO(two$ENTREZID,OrgDb = orga,pvalueCutoff = args[7],pAdjustMethod = args[8],ont = args[9])
    chebab_GO_df<-as.data.frame(chebab_GO)
    }
  
  if(args[10]=="save"){ Trim_Print(chebab_GO_df,"Enrich_result_GO") }
  
  if (args[11]=="show"){
    png("result_dot.png")
    dotplot(chebab_GO, showCategory=30) + scale_y_discrete(guide = guide_axis(n.dodge = 1))
    dev.off()
    png("result_bar.png")
    barplot(chebab_GO, showCategory=30)+scale_y_discrete(guide = guide_axis(n.dodge = 1))
    dev.off()
  }
}else{ #meaning KEGG
  filou<-read.delim2("text",header = FALSE)
  filou<-as.vector(t(filou))
  
  if(args[6]=="human"){orga="hsa"}
  else{orga="mmu"}
  
  one<-bitr_kegg(filou,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = orga)$ENSEMBL
  two<-bitr_kegg(one,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = orga)
  
  if(args[4]==1){
    fila<-read.delim2("background",header = FALSE)
    fila<-as.vector(t(fila))
    three<-bitr(fila,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = orga)$ENSEMBL
    four<-bitr(three,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = orga)
    chebab_KEGG<-enrichKEGG(two$ENTREZID,OrgDb = orga,universe =four$ENTREZID,pvalueCutoff = args[7],pAdjustMethod =args[8],ont =args[9])
    chebab_KEGG_df<-as.data.frame(chebab_KEGG)
  }
  else{
    chebab_KEGG<-enrichKEGG(two$ENTREZID,OrgDb = orga,pvalueCutoff = args[7],pAdjustMethod = args[8],ont = args[9])
    chebab_KEGG_df<-as.data.frame(chebab_KEGG)
  }
  
  if(args[10]=="save"){ Trim_Print(chebab_KEGG_df,"Enrich_result_KEGG") }
  
  if (args[11]=="show"){
    png("result_dot.png")
    dotplot(chebab_KEGG, showCategory=30) + scale_y_discrete(guide = guide_axis(n.dodge = 1))
    dev.off()
    png("result_bar.png")
    barplot(chebab_KEGG, showCategory=30)+scale_y_discrete(guide = guide_axis(n.dodge = 1))
    dev.off()
  }
}


#one=c("TIMP3","CALD1","SELENOP","MARCKS","TCF4","COL6A1","LAMC1","CST3","FSTL1","NFIB","GSN","PLAC9","LRP1","NFIA","COL1A2","TIMP2","DCN","GPX3","COL6A3")
#one<-bitr(one,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = organism)$ENSEMBL
#twooo<-bitr(one,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = organism)
#good_GO_df<-as.data.frame(enrichGO(twooo$ENTREZID,OrgDb = organism,universe = two$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",ont = "BP"))
#toTable(org.Hs.egENSEMBL2EG)$gene_id




