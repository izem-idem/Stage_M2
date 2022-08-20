
test<-all.sce.deconv[[1]]

testo<-all.sce.deconv[[2]]

test

gene_varr.test <- modelGeneVar(test)
chosen.hvgs.test <- gene_varr.test$bio > 0
gene_varr.test$chosenHvg <- chosen.hvgs.test


test <- scater::runPCA(test, subset_row = chosen.hvgs.test)
set.seed(100000)

test<-runTSNE(test,dimred="PCA")

g.test <- buildSNNGraph(test, k=10, use.dimred='PCA')
clust.test <- igraph::cluster_walktrap(g)$membership
test$label  <- factor(clust.test)

p <- plotTSNE(test, colour_by = genesToShow[3])
p <- p + ggtitle(paste("Cellules SAT", genesToShow[3]))
psat <- p

plot(psat)




lourd<-as.matrix(logcounts(current))
lourd<-as.data.frame(lourd)


write.csv(t(as.matrix(lourd[238,])),"file1.csv", row.names = TRUE)





which(rownames(test)=="CCL14")

leger<-lourd[which(rownames(test)=="APOE"),]


which(leger>0)

compare<-lourd[which(rownames(test)=="PAX7"),]

table(which(compare>0)==which(leger>0))



##investigation reconstructed
df_mnn.out<-as.data.frame(t(as.matrix(mnn.out@assays@data@listData[["reconstructed"]])))

df_mnn.out[,ncol(df_mnn.out)+1]<-mnn.out$batch

colnames(df_mnn.out)[5001]<-"Batch"

df_mnn.out<-as_tibble(t(as.matrix(lourd[238,])))

which(rownames(lourd)=="PAX7")

has_rownames(df_mnn.out)

select(df_mnn.out,c(PECAM1,Batch)) %>% filter(PECAM1 > 0.01 & Batch == 1)


######Detect how many cell sat####


current<-all.sce.deconv[[2]]

lourd<-as.matrix(logcounts(current))
lourd<-as.data.frame(lourd)

gene<-c("MYF5","PAX7")
index<-which(rownames(lourd)==gene)

df_article<-rbind(lourd[index[1],],lourd[index[2],])


df_mnn.out<-as_tibble(t(as.matrix(df_article)))


xxx<-select(df_mnn.out, c(PAX7,MYF5)) %>% filter(PAX7>0 & MYF5>0)

select(df_mnn.out, PAX7) %>% filter(PAX7>0) 

summary(xxx$"PAX7")

hist(xxx$MYF5)
