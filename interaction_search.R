#####cellphoneDB####



sce.subset=uncorrected[,c(which(uncorrected$leiden=="4"),which(uncorrected$leiden=="7"))]

sce.subset

rownames(sce.subset)<-rowData(sce.subset)$ID

maytrikss<-logcounts(sce.subset)
maytrikss<-as.data.frame(as.matrix(maytrikss))

head(maytrikss[,3])

for(i in 1:nrow(maytrikss)){
  for(j in 1:ncol(maytrikss)){
    if(maytrikss[i,j]!=0){
      maytrikss[i,j]=2^(maytrikss[i,j])
    }
  }
}

head(maytrikss[,3])

maytrikss<-cbind(rownames(maytrikss),maytrikss)

colnames(maytrikss)[1]<-"Gene"


write_tsv(maytrikss,"cpdb_counts.txt")


metadonner<-cbind(colnames(sce.subset),c(rep("Endo_cells",260),rep("Sat_cells",223)))




colnames(metadonner)<-c("Cell","cell_type")

write_tsv(as.data.frame(metadonner),"cpdb_metadata.txt")


#####CCInx#####
library(CCInx)

load("DemoDE.RData")

deL
library(readxl)

genos <- read_excel("Downloads/snRNAseq CLTI vs NI.xlsx")
genos<-as.data.frame(genos)[,1:4]


simbol<-genos$Symbol
genos<-genos[,3:4]

colnames(genos)<-c("padj","logFC")


genos$padj<-as.numeric(genos$padj)
genos$logFC<-as.numeric(genos$logFC)

row.names(genos)<-simbol


#now we have the table for sat cells let's see for the 

geny<-inspect_endo[,3:4]
colnames(geny)<-colnames(genos)
dell<-list(satcells=genos,endocells=geny)


inx<-BuildCCInx(GeneStatList=dell,
                GeneMagnitude="logFC",
                GeneStatistic="padj",
                Species="hsapiens")

PlotCCInx(INX=inx,
          cellTypeA="satcells",cellTypeB="endocells",
          proteinTypeA="Receptor",proteinTypeB="Ligand",
          TopEdges=50)

#on peut inverser had le plot ça donne 2 résulat

write_csv(inx$nodes,"inx_nodes.csv")
write_csv(inx$edges,"inx_edges.csv")


ViewCCInx(inx)


all.sce.deconv 
