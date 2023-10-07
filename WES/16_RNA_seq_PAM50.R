rm(list=ls())
options(stringsAsFactors = F)

##  setwd("../")
result_dir="16_RNA_PAM50"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}


RNA=read.csv("data/gene_tpm_anno_HXID.csv",stringsAsFactors = F,check.names = F)
RNA$gene_id=NULL


library(genefu)
data(pam50)

gene_df=pam50$centroids.map


genes_used=intersect(gene_df$EntrezGene.ID,RNA$entreID)  ## 49 genes


gene_missed=setdiff(gene_df$EntrezGene.ID,RNA$entreID)  # one gene used

gene_df[gene_df$EntrezGene.ID==gene_missed,] ## "FOXC1"


pos=match(genes_used,RNA$entreID)

PAM50_df=RNA[pos,]

rownames(PAM50_df)=PAM50_df$entreID
PAM50_df$entreID=NULL

pos=match(rownames(PAM50_df),gene_df$EntrezGene.ID)
gene_df=gene_df[pos,]

if(identical(rownames(PAM50_df),as.character(gene_df$EntrezGene.ID))){
  rownames(PAM50_df)=gene_df$probe
  
}

library(pheatmap)

pheatmap(PAM50_df,scale = "row",filename = file.path(result_dir,"01.pam50.gene.tpm.heatmap.pdf"),width = 10,height = 10)
pheatmap(log2(PAM50_df+1),scale = "row",filename = file.path(result_dir,"02.pam50.gene.log2tpm.heatmap.pdf"),width = 10,height = 10)


PAM50_df_cancer=PAM50_df[,-grep("NORM",colnames(PAM50_df))]

pheatmap(PAM50_df_cancer,scale = "row",filename = file.path(result_dir,"03.pam50.gene.cancer.tpm.heatmap.pdf"),width = 10,height = 10)
pheatmap(log2(PAM50_df_cancer+1),scale = "row",filename = file.path(result_dir,"04.pam50.gene.cancer.log2tpm.heatmap.pdf"),width = 10,height = 10)


PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(PAM50_df_cancer),
                                annot=gene_df,do.mapping=FALSE)


PAM50Preds_log<-molecular.subtyping(sbt.model = "pam50",data=t(log2(PAM50_df_cancer+1)),
                                annot=gene_df,do.mapping=FALSE)


sample_anno=data.frame(Sample=rownames(t(PAM50_df_cancer)),stringsAsFactors = F)

identical(sample_anno$Sample,names(PAM50Preds$subtype))
identical(sample_anno$Sample,names(PAM50Preds_log$subtype))



PAM50Preds_proba=as.data.frame(PAM50Preds$subtype.proba)
PAM50Preds_proba$Normal=NULL
PAM50Preds_proba$Her2=NULL


PAM50Preds_log_proba=as.data.frame(PAM50Preds_log$subtype.proba)
PAM50Preds_log_proba$Normal=NULL
PAM50Preds_log_proba$Her2=NULL



sample_anno$PAM50_TPM_subtype=apply(PAM50Preds_proba,1,function(x)names(which.max(x)))
sample_anno$PAM50_log2TPM_subtype=apply(PAM50Preds_log_proba,1,function(x)names(which.max(x)))
sample_anno$group=sapply(sample_anno$Sample,function(x)unlist(strsplit(x,"-"))[3])



sample_anno=sample_anno[order(sample_anno$PAM50_TPM_subtype),]
rownames(sample_anno)=sample_anno$Sample
sample_anno$Sample=NULL




PAM50_df_cancer_sorted=PAM50_df_cancer[,rownames(sample_anno)]

pheatmap(PAM50_df_cancer_sorted,cluster_cols = FALSE,annotation_col = sample_anno,scale = "row",
         filename = file.path(result_dir,"05.pam50.gene.cancer.tpm.heatmap.pdf"),width = 12,height = 8)

pheatmap(log2(PAM50_df_cancer_sorted+1),cluster_cols = FALSE,annotation_col = sample_anno,scale = "row"
         ,filename = file.path(result_dir,"05.pam50.gene.cancer.log2tpm.heatmap.pdf"),width = 12,height = 8)



res=t(PAM50_df_cancer_sorted)


pos=match(rownames(res),rownames(PAM50Preds_proba))

result=cbind(res,PAM50Preds_proba[pos,])


pos2=match(rownames(result),rownames(sample_anno))

result=cbind(result,sample_anno[pos2,])

result$PAM50_log2TPM_subtype=NULL

write.csv(result,file.path(result_dir,"PAM50.TPM.result.csv"),quote = F)



