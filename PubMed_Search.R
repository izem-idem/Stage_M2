library(easyPubMed)

genos=c("KLHL41", "NEB", "TCAP", "TNNT3", "TTN")

#CD9, IGFBP5, MYF5
#BTG2, FOS, MYF5, TGFBR3

#(Satellite cells) AND (myogenesis) AND (ageing) OR (senescence) AND (CD9) AND (muscle)

my_query <- '(endothelial cells) AND (HLA-C) AND (muscle)'
my_entrez_id <- get_pubmed_ids(my_query)
my_entrez_id

#fetch_all_pubmed_ids(my_entrez_id)

#my_abstracts_txt <- fetch_pubmed_data(my_entrez_id, format = "abstract",retmax = 10)
#print(my_abstracts_txt)


my_abstracts_xml <- fetch_pubmed_data(pubmed_id_list = my_entrez_id) #get them all

#my_titles <- custom_grep(my_abstracts_xml, "Abstract", "char")

#head(my_titles)
#head(my_abstracts_xml)

allFields <- table_articles_byAuth(my_abstracts_xml, included_authors = "last", getKeywords = FALSE,max_chars = 1000000)
#head(allFields)

#remove useless columns
allFields=allFields[,c(1:5)]

allFields<-as_tibble(allFields)

#remove already seen articles
#the articles that are in A and not B
my_interest<-na.omit(setdiff(allFields$pmid,like_a_deja_vu))

my_interest

allFields<-filter(allFields,pmid %in% my_interest)
allFields$year

#allFields$abstract

letter="O"

setwd("/home/izem/Desktop/STAGE/")

cat(my_query,file=paste("abstract_",letter,".txt",sep=""),sep ="\n")

#print all the abstracts so i can read them 
for(i in 1:nrow(allFields)){
  cat(paste(i,"_",i,sep=""),file =paste("abstract_",letter,".txt",sep=""),sep="\n",append =TRUE)
  cat(allFields[i,2]$doi,file = paste("abstract_",letter,".txt",sep=""),sep="\n",append =TRUE)
  cat(allFields[i,4]$abstract,file = paste("abstract_",letter,".txt",sep=""),append = TRUE,sep = "\n")  
}


#every time i read the abstract of articles i stock the pmid here so that i don't re read them another time
like_a_deja_vu<-c(like_a_deja_vu,allFields$pmid)
like_a_deja_vu<-unique(like_a_deja_vu)




