####investigations of myogenic and angiogenic markers############
angio <- read.csv("~/Downloads/angio.tsv", sep="")
angiogenesis<-read.csv("~/Downloads/angiogenesis.tsv",sep="")
angiogenic<-read.csv("~/Downloads/angiogenic.tsv",sep = "")

angiooo<-unique(rbind(angio,angiogenesis,angiogenic[,1]))

myogenic<-read.csv("~/Downloads/myogenic.tsv",sep = "")
myogenesis<-read.csv("~/Downloads/myogenesis.tsv",sep = "")

myooo<-unique(c(myogenesis[,1],myogenic[,1]))

#get index of myo and angio to get them afterward
myo_index<-which(is.element(rownames(uncorrected),myooo))
angio_index<-which(is.element(rownames(uncorrected),angiooo$Gene))


#check if the markers are among top markers genes
marker.info <- scoreMarkers(uncorrected,uncorrected$leiden)

chosen <- marker.info[["7"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
#sdd<-as.data.frame(cbind(ordered[,1:4],rownames(ordered))) # showing basic stats only, for brevity.
#ordered[1:10,1:4]



# I want two matrix (ctrl and isch) for myo and angio, and i select only satellite cells (cluster 7)
mat_angio_ctrl_sat<-as.matrix(logcounts(uncorrected)[angio_index,which(uncorrected$batch=="ctrl_2"& uncorrected$leiden=="7")])
mat_angio_isch_sat<-as.matrix(logcounts(uncorrected)[angio_index,which(uncorrected$batch=="isch"&uncorrected$leiden=="7")])

mat_myo_ctrl_sat<-as.matrix(logcounts(uncorrected)[myo_index,which(uncorrected$batch=="ctrl_2"&uncorrected$leiden=="7")])
mat_myo_isch_sat<-as.matrix(logcounts(uncorrected)[myo_index,which(uncorrected$batch=="isch"&uncorrected$leiden=="7")])

# I want two matrix (ctrl and isch) for myo and angio, and i select only endothelial cells (cluster 4 & 10)
mat_angio_ctrl_endo<-as.matrix(logcounts(uncorrected)[angio_index,which(uncorrected$batch=="ctrl_2"&uncorrected$leiden=="4" | uncorrected$leiden=="10")])
mat_angio_isch_endo<-as.matrix(logcounts(uncorrected)[angio_index,which(uncorrected$batch=="isch"&uncorrected$leiden=="4" | uncorrected$leiden=="10")])

mat_myo_ctrl_endo<-as.matrix(logcounts(uncorrected)[myo_index,which(uncorrected$batch=="ctrl_2"&uncorrected$leiden=="4" | uncorrected$leiden=="10")])
mat_myo_isch_endo<-as.matrix(logcounts(uncorrected)[myo_index,which(uncorrected$batch=="isch"&uncorrected$leiden=="4" | uncorrected$leiden=="10")])


dim(mat_myo_ctrl_sat)





#then write 4 csv
#write.csv(as.data.frame(mat_myo_ctrl),"/home/izem/Desktop/STAGE/Matrix/mat_myo_ctrl.csv",row.names = TRUE)



#get the juice baby
summary_angio_sat<-data.frame("angio_ctrl"=rowMeans(mat_angio_ctrl_sat),"angio_isch"=rowMeans(mat_angio_isch_sat),"LogFC"=rowMeans(mat_angio_ctrl_sat)-rowMeans(mat_angio_isch_sat),"percent_change"=(2^(rowMeans(mat_angio_ctrl_sat)-rowMeans(mat_angio_isch_sat)))*100)
summary_myo_sat<-data.frame("myo_ctrl"=rowMeans(mat_myo_ctrl_sat),"myo_isch"=rowMeans(mat_myo_isch_sat),"LogFC"=rowMeans(mat_myo_ctrl_sat)-rowMeans(mat_myo_isch_sat),"percent_change"=(2^(rowMeans(mat_myo_ctrl_sat)-rowMeans(mat_myo_isch_sat)))*100)

summary_angio_endo<-data.frame("myo_ctrl"=rowMeans(mat_angio_ctrl_endo),"myo_isch"=rowMeans(mat_angio_isch_endo),"LogFC"=rowMeans(mat_angio_ctrl_endo)-rowMeans(mat_angio_isch_endo),"percent_change"=(2^(rowMeans(mat_angio_ctrl_endo)-rowMeans(mat_angio_isch_endo)))*100)
summary_myo_endo<-data.frame("myo_ctrl"=rowMeans(mat_myo_ctrl_endo),"myo_isch"=rowMeans(mat_myo_isch_endo),"LogFC"=rowMeans(mat_myo_ctrl_endo)-rowMeans(mat_myo_isch_endo),"percent_change"=(2^(rowMeans(mat_myo_ctrl_endo)-rowMeans(mat_myo_isch_endo)))*100)






#####investigations of senescence markers###########

##in the  first place i will only look at p16 and p21 (CDKN1-2A)

death_index<-which(is.element(rownames(uncorrected),c("CDKN1A","CDKN2A")))


#in the second place i will look at all markers from HP atlas that have as keyword Apoptosis or Necrosis

#for Apoptosis (TSV downloaded from HP Atlas)
apop <- read.csv("~/Downloads/Apoptosis.tsv", sep="")
necro<-read.csv("~/Downloads/Necrosis.tsv", sep="")
seno<-read.csv("~/Downloads/senescence.tsv",sep="")

death_index<-which(is.element(rownames(uncorrected),seno$Gene))

mat_death_ctrl_sat<-as.matrix(logcounts(uncorrected)[death_index,which(uncorrected$batch=="ctrl_2"& uncorrected$leiden=="7")])
mat_death_isch_sat<-as.matrix(logcounts(uncorrected)[death_index,which(uncorrected$batch=="isch"&uncorrected$leiden=="7")])

mat_death_ctrl_endo<-as.matrix(logcounts(uncorrected)[death_index,which(uncorrected$batch=="ctrl_2"& uncorrected$leiden=="4" | uncorrected$leiden=="10")])
mat_death_isch_endo<-as.matrix(logcounts(uncorrected)[death_index,which(uncorrected$batch=="isch"& uncorrected$leiden=="4" | uncorrected$leiden=="10")])

#get the juice darling
summary_death_sat<-data.frame("death_ctrl"=rowMeans(mat_death_ctrl_sat),"death_isch"=rowMeans(mat_death_isch_sat),"LogFC"=rowMeans(mat_death_ctrl_sat)-rowMeans(mat_death_isch_sat),"percent_change"=(2^(rowMeans(mat_death_ctrl_sat)-rowMeans(mat_death_isch_sat)))*100)

summary_death_endo<-data.frame("death_ctrl"=rowMeans(mat_death_ctrl_endo),"death_isch"=rowMeans(mat_death_isch_endo),"LogFC"=rowMeans(mat_death_ctrl_endo)-rowMeans(mat_death_isch_endo),"percent_change"=(2^(rowMeans(mat_death_ctrl_endo)-rowMeans(mat_death_isch_endo)))*100)






###investigations of inflammatory markers########
#we'll do this investigation on Marcrophages/monoyctes cel

#i will  subset the genes that involved in inflammatory response

#get the csv. 

inflama<-read.csv("~/Downloads/inflammatory.tsv", sep="")
infla_index<-which(is.element(rownames(uncorrected),inflama$Gene))


mat_ctrl_macro<-as.matrix(logcounts(uncorrected)[infla_index,which(uncorrected$batch=="ctrl_2"& uncorrected$leiden=="2" | uncorrected$leiden=="6")])
mat_isch_macro<-as.matrix(logcounts(uncorrected)[infla_index,which(uncorrected$batch=="isch"&uncorrected$leiden=="2" | uncorrected$leiden=="6")])


summary_macro<-data.frame("Mono/Macro_ctrl"=rowMeans(mat_ctrl_macro),"Mono/Macro_isch"=rowMeans(mat_isch_macro),"LogFC"=rowMeans(mat_ctrl_macro)-rowMeans(mat_isch_macro),"percent_change"=(2^(rowMeans(mat_ctrl_macro)-rowMeans(mat_isch_macro)))*100)


######investigations of fibers###########

#let's import the genes that are involved in the secretion of muscular fiber from the article of nature "Single-cell transcriptional profiles in human skeletal muscle

library(readxl)
fyber <- read_excel("Downloads/41598_2019_57110_MOESM5_ESM.xlsx", 
                    col_types = c("text", "text", "numeric", 
                                  "numeric", "numeric"))
fyber<-fyber[-1,]

colnames(fyber)<-c("gene","muscle_type","log2 fold-change (vs. opposite fiber-type)","padj","Previous documentation as fiber-type specific")

fyber_index<-which(is.element(rownames(uncorrected),fyber$gene))

mat_ctrl_fyber<-as.matrix(logcounts(uncorrected)[fyber_index,which(uncorrected$batch=="ctrl_2"& uncorrected$leiden=="2" | uncorrected$leiden=="6")])
mat_isch_fyber<-as.matrix(logcounts(uncorrected)[fyber_index,which(uncorrected$batch=="isch"&uncorrected$leiden=="2" | uncorrected$leiden=="6")])


summary_fyber<-data.frame("Fiber_ctrl"=rowMeans(mat_ctrl_fyber),"Fiber_isch"=rowMeans(mat_isch_fyber),"LogFC"=rowMeans(mat_ctrl_fyber)-rowMeans(mat_isch_fyber),"percent_change"=(2^(rowMeans(mat_ctrl_fyber)-rowMeans(mat_isch_fyber)))*100)

summary_fyber$names<-rownames(summary_fyber)

summary_fyber<-merge(summary_fyber,fyber,by.x="names",by.y="gene")

summary_fyber<-summary_fyber[,c(-7,-8)]

myoo_ids<-bitr(myooo,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = organism)

filter(maida,gene_id %in% myoo_ids$ENTREZID,Ontology=="BP"|Ontology=="MF")

####investigate the Go::term and kegg path associated with angio markers############
tabla<-as.data.frame(toTable(org.Hs.egSYMBOL2EG))
maida<-as_tibble(as.data.frame(toTable(org.Hs.egGO)))


GOTERM[["GO:0042254"]]@Definition

#add colunm of gene names from rownames
summary_angio_sat$gene_name<-rownames(summary_angio_sat)

#convert to tibble so we can play with it
summary_angio_sat<-as_tibble(summary_angio_sat)

#we have 1 gene symbol FAM160A2
angio_ids<-bitr(summary_angio_sat$gene_name,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = organism)
lost_angio<-setdiff(summary_angio_sat$gene_name,angio_ids$SYMBOL)

#i want only the BP or MF which are associated with our gene
maida_angio<-filter(maida,gene_id %in% angio_ids$ENTREZID,Ontology=="BP"|Ontology=="MF")


#let's add the gene_name column
maida_angio<-merge(maida_angio,angio_ids,by.x="gene_id",by.y="ENTREZID")

#let's merge it with my summary_angio_sat
new_summary_angio_sat<-merge(summary_angio_sat,maida_angio,by.x="gene_name",by.y="SYMBOL")

#let's integrate the term and definition in two new columns of new_summary_angio_sat
terme[m]<-GOTERM[["GO:0042254"]]@Term
definition[m]<-GOTERM[["GO:0042254"]]@Definition

terme=NULL
definition=NULL

current<-new_summary_angio_sat #to lighten the code

#this loop takes 2 min 
for(i in 1:nrow(new_summary_angio_sat)){
  terme[i]=GOTERM[[current[i,7]]]@Term #col 7 is the column of Go_id
  definition[i]=GOTERM[[current[i,7]]]@Definition
}


new_summary_angio_sat<-cbind(current,terme,definition)
new_summary_angio_sat<-as_tibble(new_summary_angio_sat)

########Let's add KEGG also

#'arg' should be one of "Path", "Module", "ncbi-proteinid", "ncbi-geneid", "uniprot", "kegg"
kegg_keys_angio<-bitr_kegg(new_summary_angio_sat$gene_id, fromType = "ncbi-geneid", toType = "Path", organism = "hsa")
#NB : KEGG doesn't care about duplicated ids in the vector
kegg_module_angio<-bitr_kegg(new_summary_angio_sat$gene_id, fromType = "ncbi-geneid", toType = "Module", organism = "hsa")

library(KEGGREST) # if i want to use KEGG API to retrieve information
query <- keggGet("hsa04550")
query[[1]]$REL_PATHWAY


#i will just add the path ids to new_summary_angio_sat
kegg_angio_sat<-merge(distinct(new_summary_angio_sat[,c(1:6)]),kegg_keys_angio,by.x="gene_id",by.y="ncbi-geneid") #i take only col 1 to 6 to lighten the DF
kegg_angio_sat<-as_tibble(kegg_angio_sat)

query[[1]]$REL_PATHWAY


###TO print the image of the pathway
png <- keggGet("hsa04550", "image") ## retrieves the image file of a
## pathway map

setwd("~")
library(png)
library(magick)
writePNG(png, "mappa.png")
img <- magick::image_read('mappa.png')
plot(img)


#### investigate the Go::term and kegg path associated with myo markers ######

tabla<-as.data.frame(toTable(org.Hs.egSYMBOL2EG))
maida<-as_tibble(as.data.frame(toTable(org.Hs.egGO)))

#GOTERM[["GO:0042254"]]

#add colunm of gene names from rownames
summary_myo_sat$gene_name<-rownames(summary_myo_sat)

#convert to tibble so we can play with it
summary_myo_sat<-as_tibble(summary_myo_sat)

#we have 1 gene symbol FAM160A2
myo_ids<-bitr(summary_myo_sat$gene_name,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = organism)
lost_myo<-setdiff(summary_myo_sat$gene_name,myo_ids$SYMBOL)

#i want only the BP or MF which are associated with our genes
maida_myo<-filter(maida,gene_id %in% myo_ids$ENTREZID,Ontology=="BP"|Ontology=="MF")

#let's add the gene_name column
maida_myo<-merge(maida_myo,myo_ids,by.x="gene_id",by.y="ENTREZID")

#let's merge it with my summary_angio_sat
new_summary_myo_sat<-merge(summary_myo_sat,maida_myo,by.x="gene_name",by.y="SYMBOL")

#let's integrate the term and definition in two new columns of new_summary_myo_sat
terme[m]<-GOTERM[["GO:0042254"]]@Term
definition[m]<-GOTERM[["GO:0042254"]]@Definition

terme=NULL
definition=NULL

current<-new_summary_myo_sat #to lighten the code

#this loop takes 2 min 
for(i in 1:nrow(new_summary_myo_sat)){
  terme[i]=GOTERM[[current[i,7]]]@Term #col 7 is the column of Go_id
  definition[i]=GOTERM[[current[i,7]]]@Definition
}

new_summary_myo_sat<-cbind(current,terme,definition)
new_summary_myo_sat<-as_tibble(new_summary_myo_sat)

########Let's add KEGG also

#'arg' should be one of "Path", "Module", "ncbi-proteinid", "ncbi-geneid", "uniprot", "kegg"
kegg_keys_myo<-bitr_kegg(unique(new_summary_myo_sat$gene_id), fromType = "ncbi-geneid", toType = "Path", organism = "hsa")
#NB : KEGG doesn't care about duplicated ids in the vector
kegg_module_myo<-bitr_kegg(unique(new_summary_myo_sat$gene_id), fromType = "ncbi-geneid", toType = "Module", organism = "hsa")

library(KEGGREST) # if i want to use KEGG API to retrieve information
query <- keggGet("hsa04550")
query[[1]]$REL_PATHWAY


#i will just add the path ids to new_summary_angio_sat
kegg_angio_sat<-merge(distinct(new_summary_angio_sat[,c(1:6)]),kegg_keys_angio,by.x="gene_id",by.y="ncbi-geneid") #i take only col 1 to 6 to lighten the DF
kegg_angio_sat<-as_tibble(kegg_angio_sat)

query[[1]]$REL_PATHWAY


###TO print the image of the pathway
png <- keggGet("hsa04550", "image") ## retrieves the image file of a
## pathway map

setwd("~")
library(png)
library(magick)
writePNG(png, "mappa.png")
img <- magick::image_read('mappa.png')
plot(img)

#### investigate the Go::term and kegg path associated with the 30 superbes of sat cells #######

tabla<-as.data.frame(toTable(org.Hs.egSYMBOL2EG))
maida<-as_tibble(as.data.frame(toTable(org.Hs.egGO)))

#GOTERM[["GO:0042254"]]

#we'll work with ENSG ids because of some problem with Mito genes
buffer<-as_tibble(merge(buffer_index,buffer,by.x="Symbol",by.y="gene_name"))

#in order to get ids of MT- gene we should take ENSEMBL instead of gene symbol
buffer_ids<-bitr(buffer$ID,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = organism)# one ENSG have 2 ENTREZID !!! 
lost_buffer<-setdiff(buffer$ID,buffer_ids$ENSEMBL) #there is no loss of genes

#i lost one gene in translation because it simply don't have a entrezgene id !! >> AC020916.1
View(toTable(org.Hs.egENSEMBL))
#you can confirm by searching it in this database

#i want only the BP or MF which are associated with our genes
maida_buffer<-filter(maida,gene_id %in% buffer_ids$ENTREZID,Ontology=="BP"|Ontology=="MF")

#let's add the gene_name column
maida_buffer<-merge(maida_buffer,buffer_ids,by.x="gene_id",by.y="ENTREZID")

#let's merge it with my buffer
summary_buffer<-merge(buffer,maida_buffer,by.x="ID",by.y="ENSEMBL")
unique(summary_buffer$Symbol) #you can test if all the 29 genes are conserved

#let's integrate the term and definition in two new columns of buffer
terme[m]<-GOTERM[["GO:0042254"]]@Term
definition[m]<-GOTERM[["GO:0042254"]]@Definition

terme=NULL
definition=NULL

current<-summary_buffer#to lighten the code

#this loop takes 2 min 
for(i in 1:nrow(summary_buffer)){
  terme[i]=GOTERM[[current[i,9]]]@Term #col 9 is the column of Go_id
  definition[i]=GOTERM[[current[i,9]]]@Definition
}

summary_buffer<-cbind(current,terme,definition)
summary_buffer<-as_tibble(summary_buffer)


########Let's add KEGG also

#'arg' should be one of "Path", "Module", "ncbi-proteinid", "ncbi-geneid", "uniprot", "kegg"
kegg_keys_buffer<-bitr_kegg(unique(summary_buffer$gene_id), fromType = "ncbi-geneid", toType = "Path", organism = "hsa")
#NB : KEGG doesn't care about duplicated ids in the vector
kegg_module_buffer<-bitr_kegg(unique(summary_buffer$gene_id), fromType = "ncbi-geneid", toType = "Module", organism = "hsa")

library(KEGGREST) # if i want to use KEGG API to retrieve information
query <- keggGet("hsa04550")
query[[1]]$REL_PATHWAY


#i will just add the path ids to new_summary_angio_sat
kegg_buffer<-merge(distinct(summary_buffer[,1:6]),kegg_keys_angio,by.x="gene_id",by.y="ncbi-geneid") #i take only col 1 to 6 to lighten the DF
kegg_angio_sat<-as_tibble(kegg_angio_sat)

query[[1]]$REL_PATHWAY


###TO print the image of the pathway
png <- keggGet("hsa04714", "image") ## retrieves the image file of a
## pathway map

setwd("~")
library(png)
library(magick)
writePNG(png, "mappa.png")
img <- magick::image_read('mappa.png')
plot(img)



##### From GO BP search wich involved our DE genes (deprecated)############

angio_bp<-read.csv("/home/izem/Documents/angio.csv",sep = ",")
angio_bp<-as.data.frame(t(angio_bp))
real_angio_bp<-merge(angio_bp,maida,by.x="V1",by.y="go_id")


myo_bp<-read.csv("/home/izem/Documents/muscle.csv",sep = ",")
myo_bp<-as.data.frame(t(myo_bp))
real_myo_bp<-merge(myo_bp,maida,by.x="V1",by.y="go_id")



#add col Term and Def for GO
real_myo_bp$terme=rep(NA,nrow(real_myo_bp))
real_myo_bp$definition=rep(NA,nrow(real_myo_bp))

real_angio_bp$term=rep(NA,nrow(real_angio_bp))
real_angio_bp$definition=rep(NA,nrow(real_angio_bp))

#replace gene_id by gene Symbol 
#don't re run this loop it takes maybe 2 min ! 
for(i in 1:nrow(real_myo_bp)){
  real_myo_bp[i,2]=bitr(real_myo_bp[i,2],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL
  real_myo_bp[i,5]=GOTERM[[real_myo_bp[i,1]]]@Term
  real_myo_bp[i,6]=GOTERM[[real_myo_bp[i,1]]]@Definition
}

for(i in 1:nrow(real_angio_bp)){
  real_angio_bp[i,2]=bitr(real_angio_bp[i,2],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL
  real_angio_bp[i,5]=GOTERM[[real_angio_bp[i,1]]]@Term
  real_angio_bp[i,6]=GOTERM[[real_angio_bp[i,1]]]@Definition
  }

#let's add also the corresponding lFC for each gene in the 2 dfs
real_angio_bp<-merge(real_angio_bp,buffy,by.x="gene_id",by.y="gene_name")
real_myo_bp<-merge(real_myo_bp,buffy,by.x="gene_id",by.y="gene_name")

#split by up and down regulated

up_real_angio_bp<-real_angio_bp[which(real_angio_bp$summary.logFC>=0),]
down_real_angio_bp<-real_angio_bp[which(real_angio_bp$summary.logFC<0),]

up_real_myo_bp<-real_myo_bp[which(real_myo_bp$summary.logFC>=0),]
down_real_myo_bp<-real_myo_bp[which(real_myo_bp$summary.logFC<0),]


######From KEGG paths search which involve our DE genes (deprecated) ######


library(KEGGREST)
kegg_lista<-keggList("hsa")

passway <- downloadPathways("hsa")

?download_KEGG

head(sort(kegg_lista,decreasing = FALSE))

res <-keggFind("pathway", c("angiogenesis")) 




