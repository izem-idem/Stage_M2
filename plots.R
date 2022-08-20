#load Rdata and then open script_investigations_clusters.R


list_gene<-c("SPARCL1","ANXA2","MYF5","ANXA1","MEG3","ZFP36")

num_list_gene<-which(is.element(rownames(uncorrected),list_gene))

mat_myo_sat_ctrl<-as.matrix(logcounts(uncorrected)[num_list_gene,which(uncorrected$batch=="ctrl_2"& uncorrected$leiden=="7")])
mat_myo_sat_isch<-as.matrix(logcounts(uncorrected)[num_list_gene,which(uncorrected$batch=="isch" & uncorrected$leiden=="7")])

mat_myo_sat_ctrl<-t(mat_myo_sat_ctrl)
mat_myo_sat_isch<-t(mat_myo_sat_isch)

rownames(mat_myo_sat_ctrl) <- NULL
rownames(mat_myo_sat_isch) <- NULL

mat_myo_sat_ctrl<-as.data.frame(mat_myo_sat_ctrl)
mat_myo_sat_isch<-as.data.frame(mat_myo_sat_isch)


mat_myo_sat_ctrl$sample<-"ctrl"
mat_myo_sat_isch$sample<-"isch"

df_violin_myriam<-rbind(mat_myo_sat_ctrl,mat_myo_sat_isch)

mat_myo_sat_ctrl$sample<-"ctrl"

View(df_violin_myriam)

library(ggeasy)

a <- ggplot(df_violin_myriam, aes(x=sample, y=ANXA1,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("ANXA1")+
  scale_fill_brewer(palette="RdBu") + theme_classic()+
  ggeasy::easy_center_title()

b <- ggplot(df_violin_myriam, aes(x=sample, y=ANXA2,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("ANXA2")+
  scale_fill_brewer(palette="RdBu") + theme_classic()+
  ggeasy::easy_center_title()

c <- ggplot(df_violin_myriam, aes(x=sample, y=MYF5,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("MYF5")+
  scale_fill_brewer(palette="RdBu") +theme_classic()+
  ggeasy::easy_center_title()

d <- ggplot(df_violin_myriam, aes(x=sample, y=SPARCL1,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("SPARCL1")+
  scale_fill_brewer(palette="RdBu") +theme_classic()+
  ggeasy::easy_center_title()

e <- ggplot(df_violin_myriam, aes(x=sample, y=MEG3,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("MEG3")+
  scale_fill_brewer(palette="RdBu") + theme_classic()+
  ggeasy::easy_center_title()

f <- ggplot(df_violin_myriam, aes(x=sample, y=ZFP36,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("ZFP36")+
  scale_fill_brewer(palette="RdBu") + theme_classic()+
  ggeasy::easy_center_title()


gridExtra::grid.arrange(
  a,b,c,d,e,f,ncol=3
)

par(mfrow=c(1,1))

##second plot for myriam


list_gene<-c("MGP","S100A6","NME2","CST3")

num_list_gene<-which(is.element(rownames(uncorrected),list_gene))

mat_angio_sat_ctrl<-as.matrix(logcounts(uncorrected)[num_list_gene,which(uncorrected$batch=="ctrl_2"& uncorrected$leiden=="7")])
mat_angio_sat_isch<-as.matrix(logcounts(uncorrected)[num_list_gene,which(uncorrected$batch=="isch" & uncorrected$leiden=="7")])


mat_angio_sat_ctrl<-t(mat_angio_sat_ctrl)
mat_angio_sat_isch<-t(mat_angio_sat_isch)

rownames(mat_angio_sat_ctrl) <- NULL
rownames(mat_angio_sat_isch) <- NULL

mat_angio_sat_ctrl<-as.data.frame(mat_angio_sat_ctrl)
mat_angio_sat_isch<-as.data.frame(mat_angio_sat_isch)


mat_angio_sat_ctrl$sample<-"ctrl"
mat_angio_sat_isch$sample<-"isch"

df_violin_myriam<-rbind(mat_angio_sat_ctrl,mat_angio_sat_isch)

a <- ggplot(df_violin_myriam, aes(x=sample, y=MGP,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("MGP")+
  scale_fill_brewer(palette="RdBu") + theme_classic()+
  ggeasy::easy_center_title()

b <- ggplot(df_violin_myriam, aes(x=sample, y=S100A6,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("S100A6")+
  scale_fill_brewer(palette="RdBu") + theme_classic()+
  ggeasy::easy_center_title()

c <- ggplot(df_violin_myriam, aes(x=sample, y=NME2,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("NME2")+
  scale_fill_brewer(palette="RdBu") +theme_classic()+
  ggeasy::easy_center_title()

d <- ggplot(df_violin_myriam, aes(x=sample, y=CST3,fill=sample)) + 
  geom_violin()+
  xlab("")+ylab("Expression level")+ ggtitle("CST3")+
  scale_fill_brewer(palette="RdBu") +theme_classic()+
  ggeasy::easy_center_title()


gridExtra::grid.arrange(
  a,b,c,d,ncol=2
)

#######dot plot for go terms###########


#I have to add ENSEMBL data , and I also added ENTREZID data

submarine_data<-merge(subset_data_endo,as.data.frame(rowData(sce.subset_endo)),by.x="Symbol",by.y="Symbol")

submarine_data<-submarine_data[,1:6]


#lazem tu fais tourner DeSEQ.R d'abord pour avoir subset data
subhuit_data<-bitr(submarine_data$ID,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = organism)


df_go<-enrichGO(subhuit_data$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",ont = "BP",OrgDb = organism)

df_kegg<-enrichKEGG(subhuit_data$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",organism = "hsa")


View(as.data.frame(df_go))
View(as.data.frame(df_kegg))


intersect(subset_data$Symbol,subset_data_endo$Symbol)



write_csv(subset_data_endo,"DEGs_endo.csv")

View(merge(endof_dplyr,subset_data,by.x="Symbol",by.y="Symbol"))

write_csv(merge(endof_dplyr,subset_data_endo,by.x="Symbol",by.y="Symbol")[,1:2],"endo_ENSG.csv")

#Gene ratio x/y means that you have x genes over the total input gene y that are associated to this pathway
#Bg ratio A/B means that you have A genes over the total of genes B in this species that are known to be associated with this pathway.

dotplot(df_go, showCategory=20) + ggtitle("GO Biological Process") + scale_y_discrete(guide = guide_axis(n.dodge = 1))

barplot(df_kegg, showCategory=34)+scale_y_discrete(guide = guide_axis(n.dodge = 1))+ggtitle("GO Biological Process")+theme_cowplot() #scale_y_discrete to avoid overlapping in Y axis


gridExtra::grid.arrange(
  barplot(df_go, showCategory=34)+scale_y_discrete(guide = guide_axis(n.dodge = 1))+ggtitle("GO Term (BP)"), #scale_y_discrete to avoid overlapping in Y axis
  barplot(df_kegg, showCategory=34)+scale_y_discrete(guide = guide_axis(n.dodge = 1))+ggtitle("KEGG Term") #scale_y_discrete to avoid overlapping in Y axis
  ,ncol=2
)

bar_go<-barplot(df_go, showCategory=34)+scale_y_discrete(guide = guide_axis(n.dodge = 1))+ggtitle("GO Term (BP)")
bar_kegg<-barplot(df_kegg, showCategory=34)+scale_y_discrete(guide = guide_axis(n.dodge = 1))+ggtitle("KEGG Term")


bar_go+bar_kegg+plot_layout(guides = "collect")

write_csv(as.data.frame(df_go),"Tableau_GO:term_endo.csv")
write_csv(as.data.frame(df_kegg),"Tableau_KEGG:term_endo.csv")





df_go@result<-df_go@result[c(5,6,12,18,19,22,25,27,29,31,32,37,45,47,55,60,61,66,68,70,71,81,99,100,103,104,112,114,143),] # to select only the interessting terms
#te3 endo.
#c(5,6,12,18,19,22,25,27,29,31,32,37,45,47,55,60,61,66,68,70,71,81,99,100,103,104,112,114,143) 

df_kegg@result<-df_kegg@result[c(1:20,28,31,37:41),]

View(cbind(rownames(genos),genos))


lion<-merge(subset_data,subhuit_data,by.x="ENSEMBL",by.y="ENSEMBL")

for(i in 1:nrow(douf_go)){
  douf_go[i,8]=toString((bitr(strsplit(douf_go[i,8],"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL))
}


col_up=list()
col_down=list()

for(i in 1:nrow(douf_go)){ 
  vect_col_up=NULL
  vect_col_down=NULL
  x=0
  ester=str_split(douf_go[i,8],"/") #get the gene_id column and split it into a vector
  for(j in 1:length(ester[[1]])){
    x=x+1
    gene=str_trim(ester[[1]][j],"both") #remove spaces at start and end of the string
    if(lion[lion$ENTREZID==gene,4]>0){ #if gene is up regulated then 
      vect_col_up[j]=gene
    }
    if(lion[lion$ENTREZID==gene,4]<0){
      vect_col_down[j]=gene
    }
  }
  col_up[[i]]=vect_col_up
  col_down[[i]]=vect_col_down
}

GO_clone=douf_go[,c(1,2,3,6)]

GO_clone$up=rep(NA,nrow(GO_clone))
GO_clone$down=rep(NA,nrow(GO_clone))

for(i in 1:nrow(GO_clone)){
  GO_clone[i,5]=toString(na.omit(col_up[[i]]))
  #GO_clone[i,6]=toString(na.omit(col_down[[i]]))
}

GO_clone[which(is.na(GO_clone$up)),]$up=rep(GO_clone["GO:0031338",]$up,3)


for(i in 1:nrow(GO_clone)){
  GO_clone[i,5]=toString(bitr(strsplit(gsub(" ","",GO_clone[i,5]),",")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL)
  GO_clone[i,6]=toString(bitr(strsplit(gsub(" ","",GO_clone[i,6]),",")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL)
}

toString(bitr(strsplit(gsub(" ","",GO_clone[1,6]),",")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL)


write_csv(lion,"lista.csv")




########plot Cluster ################

test<-mnn.out


test$batch<-as.factor(test$batch)

levels(test$batch)<-c("CTRL","IMI")

gridExtra::grid.arrange(
  plotTSNE(test, colour_by = "Cell_types",by_exprs_values = "reconstructed"),
  plotTSNE(test, colour_by = "batch",by_exprs_values = "reconstructed"),
  ncol=2
)

levels(test$leiden)<-c("NK_cells","Monocytes","APOD_Faps","Endothelial_cells","Faps","Macrophages","MuSCs","CMP","T_cells","PCV_Endothelial_cells","Pericytes")

colnames(colData(test))[15]<-"Cell_types"



