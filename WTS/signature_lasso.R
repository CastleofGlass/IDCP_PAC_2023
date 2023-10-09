

setwd("signature/")

library(stringr)
library(dplyr)
library(pROC)
library(openxlsx)
library(glmnet)

#pair==============================================================
pairsample <- c('00.IDCP','01.IDCP','02.IDCP','03.IDCP','04.IDCP','05.IDCP','06.IDCP','07.IDCP','08.IDCP','09.IDCP','10.IDCP','11.IDCP','13.IDCP','14.IDCP','17.IDCP',
                '00.PCA','01.PCA','02.PCA','03.PCA','04.PCA','05.PCA','06.PCA','07.PCA','08.PCA','09.PCA','10.PCA','11.PCA','13.PCA','14.PCA','17.PCA')
deseq_pair <- read.xlsx("IDCP_vs_PCA_pair差异分析汇总.xlsx",sheet = 1)
diffGeneListU <- deseq_pair[which(deseq_pair$pvalue<0.01 & deseq_pair$log2FoldChange > 2),"gene"] %>% unique()
diffGeneListD <- deseq_pair[which(deseq_pair$pvalue<0.01 & deseq_pair$log2FoldChange < -2),"gene"] %>% unique()
diffGeneList <- c(diffGeneListU,diffGeneListD)
#nopair============================================================
# nopairsample <- 
deseq <- read.table("IDCP_vs_PCA_deseq2_gene.tsv",sep = "\t",header = T,stringsAsFactors = F)
diffGeneListU <- deseq[which(deseq$pvalue<0.01 & deseq$log2FoldChange > 2),] %>% rownames()
diffGeneListD <- deseq[which(deseq$pvalue<0.01 & deseq$log2FoldChange < -2),] %>% rownames()
diffGeneList <- c(diffGeneListU,diffGeneListD)
# diffGeneList <- read.table("diffList.txt",sep = "\t",header = F,stringsAsFactors = F)

#cell cycle 20211030=======================================================
#go bp
go_bp <- read.xlsx("IDCP_vs_PCA_pair差异分析汇总.xlsx",sheet = 2)
go_bp_list <- go_bp[grep("cell cycle",go_bp$Description,ignore.case = T),"core_enrichment"] %>% str_split("\\/") %>% unlist() %>% unique()
kegg <- read.xlsx("IDCP_vs_PCA_pair差异分析汇总.xlsx",sheet = 4)
kegg_list <- kegg[grep("cell cycle",kegg$Description,ignore.case = T),"core_enrichment"] %>% str_split("\\/") %>% unlist() %>% unique()

go_bp <- read.xlsx("IDCP_vs_PCA差异分析汇总.xlsx",sheet = 2)
go_bp_list <- go_bp[grep("cell cycle",go_bp$Description,ignore.case = T),"core_enrichment"] %>% str_split("\\/") %>% unlist() %>% unique()
kegg <- read.xlsx("IDCP_vs_PCA差异分析汇总.xlsx",sheet = 4)
kegg_list <- kegg[grep("cell cycle",kegg$Description,ignore.case = T),"core_enrichment"] %>% str_split("\\/") %>% unlist() %>% unique()


cellcycle_gene <- unique(c(go_bp_list,kegg_list))
library(org.Hs.eg.db)
library(clusterProfiler)
# 
IDmapping <- bitr(cellcycle_gene,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db)
# gene_info <- read.table("gene_infor.xls",sep = "\t",stringsAsFactors = F)
IDmapping$ENSEMBL %in% str_remove(geneTPM$gene_id,"\\..*")

diffGeneList <- c(str_remove(diffGeneList,"\\..*"),IDmapping$ENSEMBL) %>% unique()

geneTPM <- read.table(file = "gene_tpm_anno_HXID.txt",sep = "\t",header = T,stringsAsFactors = F)
# geneTPM <- geneTPM[apply(geneTPM[,2:50],1,function(x) ifelse(max(x)>=2,TRUE,FALSE)),]
colnames(geneTPM) <- str_remove(colnames(geneTPM),"HXIDCP\\.")

# rownames(geneTPM) <- geneTPM$gene_id
geneTPM$gene_id <- str_remove(geneTPM$gene_id,"\\..*")

diffTPM <- geneTPM[geneTPM$gene_id %in% diffGeneList,
                         c(1,grep("IDCP",colnames(geneTPM)),grep("PCA",colnames(geneTPM)))] %>% unique()

#paired
diffTPM <- geneTPM[geneTPM$gene_id %in% diffGeneList,
                   c("gene_id",pairsample)] %>% unique()


rownames(diffTPM) <- diffTPM$gene_id

#-----------------------------------------------


#tpm在每一组中的log均值需要大于1或2或4，实验才能验证

# diffTPM <- diffTPM[which(apply(diffTPM[,grep("IDCP",colnames(diffTPM))],1,mean)>1 & apply(diffTPM[,grep("PCA",colnames(diffTPM))],1,mean)>1),]
# diffTPM <- diffTPM[which(apply(diffTPM[,grep("IDCP",colnames(diffTPM))],1,mean)>=2 | apply(diffTPM[,grep("PCA",colnames(diffTPM))],1,mean)>=2),]

mean1_idcp <- apply(diffTPM[,grep("IDCP",colnames(diffTPM))],1,mean)
mean1_pca <- apply(diffTPM[,grep("PCA",colnames(diffTPM))],1,mean)

diffTPM <- diffTPM[which(mean1_idcp>1 & mean1_pca>1),]

#tpm在每一组中非0值样本超过一半
num0_idcp <- apply(diffTPM[,grep("IDCP",colnames(diffTPM))],1,function(x){sum(x > 0)})
num0_pca <- apply(diffTPM[,grep("PCA",colnames(diffTPM))],1,function(x){sum(x > 0)})

diffTPM <- diffTPM[which(num0_idcp>(sum(grepl("IDCP",colnames(diffTPM)))/2) & num0_pca>(sum(grepl("PCA",colnames(diffTPM)))/2)),]

#tpm在每一组中cv小于
cv_idcp <- apply(diffTPM[,grep("IDCP",colnames(diffTPM))],1,function(x){x1 <- x[x>0]; sd(x1)/mean(x1)})
cv_pca <- apply(diffTPM[,grep("PCA",colnames(diffTPM))],1,function(x){x1 <- x[x>0]; sd(x1)/mean(x1)})

diffTPM <- diffTPM[which(cv_idcp < 0.4 & cv_pca < 0.4),]


#cv值--------------------------------
num_cv <- function(x){
  which(cv_idcp < x & cv_pca < x) %>% length()
}

plot(seq(0,3,0.01),sapply(seq(0,3,0.01),num_cv),xlab="CV",ylab="Number of gene",pch=16,col=rgb(.1,.7,0.6,alpha=0.6))

#------------------------------------

exp_lasso <- as.data.frame(t(diffTPM[,-1]))
marker <- diffTPM$gene_id
exp_lasso <- as.data.frame(t(diffTPM[genelogi,-1]))

# exp_lasso$Sample <- rownames(exp_lasso)
exp_lasso$type <- str_remove(rownames(exp_lasso),"\\d+\\.")
exp_lasso$type[exp_lasso$type=="IDCP"] <- 1
exp_lasso$type[exp_lasso$type=="PCA"] <- 0
exp_lasso$type <- as.numeric(exp_lasso$type)
lx <- model.matrix(object = type ~ ., data = exp_lasso)[, -1]
ly <- exp_lasso$type

cv_lasso <- cv.glmnet(x = lx, y = ly, family="gaussian",
          nlambda=200,alpha = 1)
plot(cv_lasso)

l_pred <- predict(cv_lasso,newx = lx, s = cv_lasso$lambda.min)

pred_res <- as.numeric(l_pred > median(l_pred))

exp_lasso_pred_res <- exp_lasso
exp_lasso_pred_res$pred <- l_pred

lroc <-roc(exp_lasso_pred_res$type, exp_lasso_pred_res$pred)

# table(Batch2_y);coef(Batch2_cv_lasso, s=Batch2_cv_lasso$lambda.min)
table(ly)
cv_lasso$lambda.min
tmp <- coef(cv_lasso, s=cv_lasso$lambda.min)
tmp
diffTPM[genelogi,]$gene_id[tmp@i]

plot(lroc,print.auc=TRUE,plot=TRUE,print.thres=TRUE)

# gene_info <- read.table("../gene_infor.xls",sep = "\t",header = T,stringsAsFactors = F)
# 
# gene_info <- filter(gene_info,entreID!="no")
# 
# sum(unique(gene_info$gene_id) %in% diffGeneList$V1)
# 
# sum(str_remove(unique(gene_info$gene_id),"\\..*") %in% str_remove(diffGeneList$V1,"\\..*"))


diffTPM

diffTPM_anno <- left_join(unique(diffTPM),unique(gene_info[,3:7]),by=c("gene_id"="gene_id"))
# diffTPM_anno <- filter(diffTPM_anno,!is.na(gene_name))
IDmapping2 <- bitr(diffTPM[genelogi,]$gene_id[tmp@i],fromType = "ENSEMBL",toType = c("ENTREZID","SYMBOL"),OrgDb = org.Hs.eg.db)

kegg_res <- enrichKEGG(c(IDmapping2$ENTREZID,107075247),organism = "hsa",keyType = "kegg",pvalueCutoff = 1,pAdjustMethod = "BH")
go_bp_res <- enrichGO(c(IDmapping2$ENTREZID,107075247),keyType = "ENTREZID",OrgDb = 'org.Hs.eg.db',ont = "BP",pvalueCutoff = 1,pAdjustMethod = "BH")
# gene_info2 <- read.table("../gene2ensembl_210902.qinc",sep = "\t",header = T,stringsAsFactors = F)
# sum(unique(gene_info2$ensembl) %in% str_remove(diffGeneList$V1,"\\..*"))

barplot(go_bp_res,showCategory=20,title="Biological Process",)
barplot(kegg_res,showCategory=20,title="KEGG",)

write.table(go_bp_res@result,"logistic_res_gobp.txt",sep = "\t",row.names = F,quote = F)
write.table(kegg_res@result,"logistic_res_kegg.txt",sep = "\t",row.names = F,quote = F)

# CV<1
# ENSG00000080839  0.028307056
# ENSG00000112312  0.013602632
# ENSG00000113300  0.002452122
# ENSG00000134480 -0.003036497
# ENSG00000139351  0.011193525
# ENSG00000187790  0.064700742

# CV
# ENSG00000080839  0.008553443
# ENSG00000095002  0.003731308
# ENSG00000112312  0.005828633
# ENSG00000187790  0.017579378
# ENSG00000206172 -0.002189042

# CV<0.35
# ENSG00000100393  0.02297748
# ENSG00000135679  0.01208971
# ENSG00000177302  0.03717709

plot(cv_idcp,deseq[str_remove(rownames(deseq),"\\..*") %in% names(cv_idcp),"log2FoldChange"],pch=16,col=rgb(.1,.7,0.6,alpha=0.6),xlab="CV IDCP&PCA",ylab="log2FC")
points(cv_pca,deseq[str_remove(rownames(deseq),"\\..*") %in% names(cv_pca),"log2FoldChange"],pch=16,col=rgb(.6,.2,0.6,alpha=0.6))

points(cv_idcp[intersect(names(cv_idcp),IDmapping$ENSEMBL)],deseq[str_remove(rownames(deseq),"\\..*") %in% intersect(names(cv_idcp),IDmapping$ENSEMBL),"log2FoldChange"],pch=17,col="green")
points(cv_pca[intersect(names(cv_pca),IDmapping$ENSEMBL)],deseq[str_remove(rownames(deseq),"\\..*") %in% intersect(names(cv_pca),IDmapping$ENSEMBL),"log2FoldChange"],pch=17,col="red")

#
plot(cv_idcp,deseq_pair[str_remove(deseq_pair$gene,"\\..*") %in% names(cv_idcp),"log2FoldChange"],pch=16,col=rgb(.1,.7,0.6,alpha=0.6),xlab="CV IDCP&PCA",ylab="log2FC")
points(cv_pca,deseq_pair[str_remove(deseq_pair$gene,"\\..*") %in% names(cv_pca),"log2FoldChange"],pch=16,col=rgb(.6,.2,0.6,alpha=0.6))

points(cv_idcp[intersect(names(cv_idcp),IDmapping$ENSEMBL)],deseq_pair[str_remove((deseq_pair$gene),"\\..*") %in% intersect(names(cv_idcp),IDmapping$ENSEMBL),"log2FoldChange"],pch=17,col="green")
points(cv_pca[intersect(names(cv_pca),IDmapping$ENSEMBL)],deseq_pair[str_remove((deseq_pair$gene),"\\..*") %in% intersect(names(cv_pca),IDmapping$ENSEMBL),"log2FoldChange"],pch=17,col="red")


plot(cv_idcp,cv_pca)
##########
