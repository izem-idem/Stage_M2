library(org.Hs.eg.db)
library(GO.db)


#buffer_up<-inspect_sat[inspect_sat$summary.logFC>0.9,]
#buffer_down<-inspect_sat[inspect_sat$summary.logFC<=-0.9,]

#get the genes that are above 0.9 and under -0.9
buffer<-inspect_sat[inspect_sat$summary.logFC>0.9 | inspect_sat$summary.logFC<=-0.9,]

buffer<-as_tibble(buffer)


#since i saw that you should have a lot of genes and they don't have to be with good pvalue 
buffy<-inspect_sat[inspect_sat$summary.logFC!=0,]

nrow(inspect_sat[inspect_sat$summary.logFC>=0.5 | inspect_sat$summary.logFC<=-0.5,])

#i had some problems so i have to modify some of them so they can match i go check the ensg id in the web and then i get the gene_id
current<-uncorrected #take whatever sce

#from GENESYMBOL to ENSG num (no need)
buffer_index<-rowData(current)[which(is.element(rowData(current)[[2]],buffer$gene_name)),c(1,2)]

buffer_index<-as_tibble(buffer_index)


#i don't know why but i have a 10 genes which are lost in the conversions
buffy_index<-select(filter(as_tibble(rowData(current)),Symbol %in% buffy$gene_name),ID,Symbol)

current_tibble<-as_tibble(rowData(current))

#check which are these missings genes
lost_in_translate<-which(!is.element(buffy$gene_name,buffy_index$Symbol))

lost_in_translate<-buffy[lost_in_translate,]$gene_name
#i found that 11 genes have this form  TBCE_ENSG00000284770 and so to keep only the symbol do gsub("_ENSG[0-9]*","","TBCE_ENSG00000284770")
filter(buffy,gene_name %in% lost_in_translate )
#i saw in filter that they don't have much relevance so let them get lost

#select(filter(as_tibble(rowData(current)),Symbol %in% buffer$gene_name),ID)==buffer_index
#to check if buffer_index is ok


#to find gene by ENS
nomgene <- toTable(org.Hs.egENSEMBL)

#buffer_nom<-nomgene$gene_id[which(is.element(nomgene$ensembl_id,buffer_index))]
#ENSG00000269028 is present two times so it has 2 gene_id

nomgene<-as_tibble(nomgene)

tabla<-filter(nomgene,ensembl_id %in% buffer_index)

maida<-as.data.frame(toTable(org.Hs.egGO))

GOGOGO<-merge(tabla, maida, by.x="gene_id", by.y="gene_id")


#stock it into 3 different dataframes
DEgenes_MF<-GOGOGO[which(GOGOGO$Ontology=="MF"),]
DEgenes_CC<-GOGOGO[which(GOGOGO$Ontology=="CC"),]
DEgenes_BP<-GOGOGO[which(GOGOGO$Ontology=="BP"),]

gros_bucket<-list(DEgenes_BP,DEgenes_CC,DEgenes_MF)

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


DEgenes_BP<-gros_bucket[[1]]
DEgenes_CC<-gros_bucket[[2]]
DEgenes_MF<-gros_bucket[[3]]

#add the gene symbol

DEgenes_BP<-merge(select(filter(as_tibble(rowData(current)),Symbol %in% buffer$gene_name),ID,Symbol),DEgenes_BP,by.x="ID", by.y="ensembl_id")
DEgenes_CC<-merge(select(filter(as_tibble(rowData(current)),Symbol %in% buffer$gene_name),ID,Symbol),DEgenes_CC,by.x="ID", by.y="ensembl_id")
DEgenes_MF<-merge(select(filter(as_tibble(rowData(current)),Symbol %in% buffer$gene_name),ID,Symbol),DEgenes_MF,by.x="ID", by.y="ensembl_id")

#add information about logFC

DEgenes_BP<-merge(DEgenes_BP,inspect_sat,by.x="Symbol",by.y="gene_name")
DEgenes_CC<-merge(DEgenes_CC,inspect_sat[,-5],by.x="Symbol",by.y="gene_name")
DEgenes_MF<-merge(DEgenes_MF,inspect_sat[,-5],by.x="Symbol",by.y="gene_name")


library(icesTAF)
mkdir("/home/izem/Desktop/STAGE/DE_genes_GO")

setwd("/home/izem/Desktop/STAGE/")

write_csv(DEgenes_BP,paste("BP_DEgenes",".csv",sep = ""))
write_csv(DEgenes_CC,paste("CC_DEgenes",".csv",sep = ""))
write_csv(DEgenes_MF,paste("MF_DEgenes",".csv",sep = ""))


write_csv(endof_dplyr,"Touslesgénes.csv")

write_csv(GO_df,"Enrichissement_Fonctionnel_GO.csv")

write_csv(circ_input_two,"ok.csv")





#### TopGO explore (it's all in Biostar)########
library(topGO)
library(KEGGREST)

require(org.Hs.eg.db)


buffy_real_index<-select(filter(nomgene,ensembl_id %in% as.vector(buffy_index$ID)),gene_id,ensembl_id)
#AGAIN !! i had a loss of information during the translation it better not be important genes

lost_again<-which(!is.element(buffy_index$ID,buffy_real_index$ensembl_id))
lost_again<-buffy_index[lost_again,]

View(filter(buffy,gene_name %in% lost_again$Symbol))
#this genes are not present in org.Hs.egENSEMBL tout simplement


tablita<-as_tibble(toTable(org.Hs.egSYMBOL2EG))

#I select the real gene_symbol used by GO ontology and then used them for Go enrichment analysis
buffer_real_index<-select(filter(tablita,gene_id %in% tabla$gene_id[-31]),symbol)



genes<- buffy$p.value    #buffer$p.value
names(genes)<-buffy$gene_name          #buffer_real_index$symbol


allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=10)



#The p-values are used to rank the genes, which is important when using the Kolmogorov-Smirnov test.
results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")

#"Annotated" is the number of genes in your database which are annotated with the GO-term; 
#"Significant" is the number of genes from your input which are annotated with the GO-term. 
#"Expected" is an estimate of the number of genes a node of size Annotated would have if the significant genes were to be randomly selected from the gene universe.
goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=50)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.01,] #or 0.05
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
#since the . has a special meaning kevin blighe wanted to tell i really want to grab an actual point not "any characters"
#(meaning of the . in regex) , he does it the first time but he wanted to be sure to rm all the ... in a string
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

require(ggplot2)
ggplot(goEnrichment, aes(x=Term, y=-log10(KS),fill=Term)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle("Title") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw(base_size=5) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=5, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1,"cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=10),  #Text size
    title=element_text(size=15)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

#######Dot-plot#####
library(scales)

ntop <- 50
ggdata <- goEnrichment[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
gg1 <- ggplot(ggdata,
              aes(x = Term, y = -log10(KS), size = -log10(KS), fill = -log10(KS))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO Biological processes',
    subtitle = paste('Top', ntop ,'terms ordered by Kolmogorov-Smirnov p-value'),
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             size = c(0.5, 1.5, 3)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

#######KEGG pathways##############

library(clusterProfiler)

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

#stop playing with dplyr

endof_dplyr<-merge(dedup_ids,buffy_index,by.x="ENSEMBL",by.y="ID")

#since buffy and dedup_ids do not have the same length i'll have to synchronize buffy in order to have the same length
eof_dplyr<-filter(buffy,gene_name %in% endof_dplyr$Symbol)

#now i defenetly stop playing with dplyr
endof_dplyr<-merge(endof_dplyr,eof_dplyr,by.x="Symbol",by.y="gene_name")


keytypes(org.Hs.eg.db)


Go_gene_list<-endof_dplyr$summary.logFC

names(Go_gene_list)<-endof_dplyr$ENSEMBL

Go_gene_list = sort(Go_gene_list, decreasing = TRUE)



#ont one of “BP”, “MF”, “CC” or “ALL”
#nPerm the higher the number of permutations you set, the more accurate your result will, but the longer the analysis will take.
#minGSSize minimum number of genes in set (gene sets with lower than this many genes in your dataset will be ignored).
#maxGSSize maximum number of genes in set (gene sets with greater than this many genes in your dataset will be ignored).
#pvalueCutoff pvalue Cutoff.
#pAdjustMethod one of “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”

organism = "org.Hs.eg.db"


###KEGG not yet (repeat GO)
gse <- gseGO(geneList=Go_gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism)

View(as.data.frame(gse@result))

#check if genes are in the KEGG pathway or KEGG module 
#bitr_kegg(tabla$gene_id, fromType = "kegg", toType = "Path", organism = "hsa")
#bitr_kegg(tabla$gene_id, fromType = "kegg", toType = "Module", organism = "hsa")


require("DOSE")
dotplot(gse, showCategory=10, split=".sign")


View(as.data.frame(gse))

# Create a vector with the LogFC
kegg_gene_list <- endof_dplyr$summary.logFC


# Name vector with ENTREZ ids
names(kegg_gene_list) <- endof_dplyr$ENTREZID

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

# Gene ratio (count of core enrichment genes) / (count of pathway genes)  

View(as.data.frame(kk2))


library(pathview)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa02010", species = kegg_organism,limit = list(gene=5, cpd=1))

knitr::include_graphics("hsa02010.pathview.png")


####Method enrichKEGG (also the one studied in M2, i mean the ORA )########

buffalo<-inspect_sat[inspect_sat$summary.logFC>=0.5 | inspect_sat$summary.logFC<=-0.5,]

kind_buffalo<-inspect_sat[inspect_sat$summary.logFC>=0.4 | inspect_sat$summary.logFC<=-0.4,]

#Let's consider DE genes as genes having > 0.5 or < 0.5 LFC 

current<-uncorrected 

#from GENESYMBOL to ENSG num (no need) , i know i could've use endof_dplyr to get the entrezID
ENSG_buffalo<-as.data.frame(rowData(current)[which(is.element(rowData(current)[[2]],buffalo$gene_name)),c(1,2)])
buffalo<-merge(ENSG_buffalo,buffalo,by.x="Symbol",by.y="gene_name")

ENSG_kind_buffalo<-as.data.frame(rowData(current)[which(is.element(rowData(current)[[2]],kind_buffalo$gene_name)),c(1,2)])
kind_buffalo<-merge(ENSG_kind_buffalo,kind_buffalo,by.x="Symbol",by.y="gene_name")


#remove one column 
buffalo<-buffalo[,-7]
kind_buffalo<-kind_buffalo[,-7]

#now let's pass to from ENSEMBL to ENTREZID (lost 1 in translation)

buffalo_ids<-bitr(buffalo$ID,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = organism)

kind_buffalo_ids<-bitr(kind_buffalo$ID,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = organism)
  
#add col ENTREZID to the buffalo
buffalo<-merge(buffalo,buffalo_ids,by.x="ID",by.y="ENSEMBL")

kind_buffalo<-merge(kind_buffalo,kind_buffalo_ids,by.x="ID",by.y="ENSEMBL")
  
#View(as.data.frame(enrichKEGG(kind_buffalo_ids$ENTREZID,organism = kegg_organism,universe = dedup_ids$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH")))

wealthy_kegg<-enrichKEGG(buffalo_ids$ENTREZID,organism = kegg_organism,universe = dedup_ids$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH")
kegg_df<-as.data.frame(wealthy_kegg)

dotplot(wealthy_kegg, showCategory=20,split=".sign") + facet_grid(.~.sign)

barplot(wealthy_kegg, showCategory=15)

#convert count column with gene symbol
for(i in 1:nrow(kegg_df)){
  kegg_df[i,8]=toString((bitr(strsplit(kegg_df[i,8],"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL))
}


#use this library to help visualize your result
library(googleVis)
gvt = gvisTable(kegg_df)
plot(gvt)

#let's see on which pathways ORA et GSEA agrees upon
agree<-intersect(kk2@result$Description,kegg_df$Description)
#the pathways that are in A and not B 
setdiff(kk2@result$Description,kegg_df$Description)
setdiff(kegg_df$Description,kk2@result$Description)

####Method enrichGO (also the one studied in M2, i mean the ORA ) for sat cells########

kind_GO_df<-as.data.frame(enrichGO(kind_buffalo_ids$ENTREZID,OrgDb = organism,universe = dedup_ids$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",ont = "BP"))

?enrichGO
wealthy_GO<-enrichGO(buffalo_ids$ENTREZID,OrgDb = organism,universe = dedup_ids$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",ont = "BP")
GO_df<-as.data.frame(wealthy_GO)

#little warning : i gave 85 input genes put in the gene ratio they show only 83

#Gene ratio x/y means that you have x genes over the total input gene y that are associated to this pathway
#Bg ratio A/B means that you have A genes over the total of genes B in this species that are known to be associated with this pathway.

dotplot(wealthy_GO, showCategory=20) + ggtitle("GO Biological Process") + scale_y_discrete(guide = guide_axis(n.dodge = 1))
barplot(wealthy_GO, showCategory=20)+scale_y_discrete(guide = guide_axis(n.dodge = 1))+ggtitle("GO Biological Process") #scale_y_discrete to avoid overlapping in Y axis

for(i in 1:nrow(GO_df)){
  GO_df[i,8]=toString((bitr(strsplit(GO_df[i,8],"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL))
}

for(i in 1:nrow(kind_GO_df)){
  kind_GO_df[i,8]=toString((bitr(strsplit(kind_GO_df[i,8],"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL))
}


#another to visualise the result
library(GOplot)
#see the example in help( ?circle_dat) to get what kind of input you should put

circ_input_one<-cbind(GO_df[,c(1,2,6,8,9)],rep("BP",nrow(GO_df)))
#we can also do the chord chart for kegg pathways 
cirkegg_input_one<-cbind(kegg_df[,c(1,2,6,8,9)],rep("BP",nrow(kegg_df)))
colnames(circ_input_one)[6]<-c("Category")
#for kegg
colnames(cirkegg_input_one)[6]<-c("Category")

#reorder the columns
circ_input_one<-circ_input_one[,c(6,1,2,4,3,5)]
colnames(circ_input_one)<-colnames(EC$david)
#rename the 6th col to "count"
colnames(circ_input_one)[6]<-"count"

#reorder the columns for kegg
cirkegg_input_one<-cirkegg_input_one[,c(6,1,2,4,3,5)]
colnames(cirkegg_input_one)<-colnames(EC$david)
#rename the 6th col to "count"
colnames(cirkegg_input_one)[6]<-"count"




circ_input_two<-buffalo[,c(7,6)]
#from entrez_id to symbol
for(i in 1:nrow(circ_input_two)){
  circ_input_two[i,1]=bitr(circ_input_two[i,1],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL
}
colnames(circ_input_two)<-c("ID","logFC")

#inverse the value so that marianne can understand
circ_input_two$logFC<-(circ_input_two$logFC)*-1


circ<-circle_dat(circ_input_one,circ_input_two)
cirkegg<-circle_dat(cirkegg_input_one,circ_input_two)

chord<-chord_dat(data = circ,genes = circ_input_two, process = circ_input_one$Term[muscle_index[21:length(muscle_index)]]) #don't choose all the process 
#For kegg
chord_kegg<-chord_dat(data = cirkegg,genes = circ_input_two, process = cirkegg_input_one$Term)

#here i can choose to highlight only the process associated with muscle
#library(randomcoloR) #to get random colors
GOChord(chord_kegg,space = 0.01,gene.order = "logFC",
        limit = c(0,0),
        gene.space = 0.25,
        gene.size = 5,
        border.size = 0.1,
        process.label = 10)


#put num so i can get them
rownames(circ_input_one)<-seq(1,nrow(circ_input_one))

muscle_index<-c(22,27,29,31,32,33,35,38,45,47,50,60,62,63,64,66,72,89,91,93,95,102,103,109,110,67)
blood_et_vessel_index<-c(36,37,39,40,41,43,46,52,53,54,56,57,58,59,61)
energy_index<-c(1:9,13,16,24,71,23,70,10)
ionic_index<-c(11,12,14,15,19,20,21,25,26,28,30,34,42,49)
response_index<-c(17,44,51,55,65,69,73:75,81:83,86,90)
Nar_index<-c(107,105,104,78,79)

#I curated 90 BP only 20 remains
View(circ_input_one[c(-muscle_index,-blood_et_vessel_index,-energy_index,-ionic_index,-response_index,-Nar_index),])

View(circ_input_one[c(muscle_index),])


#another way to visualise.

l1 <- subset(circ, term == 'myotube differentiation', c(genes,logFC))
l2 <- subset(circ, term == 'oxidative phosphorylation', c(genes,logFC))
l3 <- subset(circ, term == 'response to copper ion', c(genes,logFC))

GOVenn(l1,l2,l3,label = c('oxidative phosphorylation','myotube differentiation'))


?GOVenn
#another way to visualise the results
parwise_GO<-enrichplot::pairwise_termsim(wealthy_GO)
emapplot(parwise_GO)



#Another way to visualise the results
cnetplot(parwise_GO)

#same as wealthy_go but SYmbol instead of EntrezID
cnet_GO<-setReadable(wealthy_GO,'org.Hs.eg.db', 'ENTREZID')

input_FC_cnet<-circ_input_two$logFC
names(input_FC_cnet)<-circ_input_two$ID
cnetplot(cnet_GO, foldChange=input_FC_cnet, circular = TRUE, colorEdge = TRUE,showCategory = 10) 

View(cnet_GO@result)
#another way to visualise your dataframe
gvt = gvisTable(GO_df)
plot(gvt)


###enrich DAVID########

library()
wealthy_DAVID<-enrichDAVID(buffalo_ids$ENTREZID,
                           universe = dedup_ids$ENTREZID,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           annotation = "GOTERM_BP_FAT")
#RDAVIDWebService is deprecated so we don't have enrich DAVID




###Create the data for enrich GO/KEGG with endothelial cells#######

lion<-inspect_endo[inspect_endo$summary.logFC>=0.4 | inspect_endo$summary.logFC<=-0.4,]

filtro<-filter(endof_dplyr,Symbol %in% lion$gene_name)

filtro<-filtro[,c(1:3)]

forgotten_genes<-setdiff(lion$gene_name,filtro$Symbol)
#let's remove HELLPAR
forgotten_genes<-forgotten_genes[-3]

#the ensembl id of the forgotten genes
remembered_genes<-bitr(c("ENSG00000141469","ENSG00000237433","ENSG00000236144","ENSG00000125810"),fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = organism)
remembered_genes<-cbind(remembered_genes,forgotten_genes)

subterfuge<-cbind(filter(lion,gene_name %in% remembered_genes$forgotten_genes),remembered_genes)

lion<-merge(filtro,inspect_endo,by.x="Symbol",by.y="gene_name")


colnames(subterfuge)
colnames(lion)

subterfuge<-subterfuge[,c(-9)]
subterfuge<-subterfuge[,c(6,7,8,1:5)]

colnames(subterfuge)<-colnames(lion)

#finally
lion<-rbind(lion,subterfuge)

#remove the gene that have FDR higher than 0.05
lion<-filter(lion,FDR<0.05)

####Method enrichGO/KEGG (ORA ) for endothelial cells##########
endo_GO_df<-as.data.frame(enrichGO(lion$ENTREZID,OrgDb = organism,universe = dedup_ids$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",ont = "BP"))


#Gene ratio x/y means that you have x genes over the total input gene y that are associated to this pathway
#Bg ratio A/B means that you have A genes over the total of genes B in this species that are known to be associated with this pathway.

dotplot(endo_GO_df, showCategory=30) + ggtitle("GO Biological Process")
barplot(endo_GO_df, showCategory=50)+scale_y_discrete(guide = guide_axis(n.dodge = 1)) #scale_y_discrete to avoid overlapping in Y axis

#change gene ids into gene names
for(i in 1:nrow(endo_GO_df)){
  endo_GO_df[i,8]=toString((bitr(strsplit(endo_GO_df[i,8],"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism)$SYMBOL))
}

## i have deleted all the tools of visualization but i can easily find them in the enrichGO section for sat cells above !
#another way to visualise your dataframe
gvt = gvisTable(GO_df)
plot(gvt)

###How about kegg pathway

endo_kegg_df<-as.data.frame(enrichKEGG(lion$ENTREZID,organism = kegg_organism,universe = dedup_ids$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH"))

######Method GSEA for endothelial cells#####


for_gene_list<-merge(as.data.frame(rowData(uncorrected))[,c(1,2)],inspect_endo,by.x="Symbol",by.y="gene_name")
for_gene_list<-merge(bitr(for_gene_list$ID,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = organism),for_gene_list,by.x="ENSEMBL",by.y="ID")

endo_gene_list<-(for_gene_list$summary.logFC)*-1 #in order to put up regulated gene at the top of the ranked list
names(endo_gene_list)<-for_gene_list$ENTREZID
endo_gene_list<-sort(endo_gene_list,decreasing = TRUE)


gse <- gseGO(geneList=endo_gene_list, 
             ont ="BP", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism,pAdjustMethod = "BH")

endo_gene_list

View(gse@result)

gse_endo <- gseKEGG(geneList= endo_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")

View(gse_endo@result)
?gseKEGG
?gseGO
