#  setwd("../")


rm(list=ls())

outdir="12_gistic_broad"

if(!dir.exists(outdir)){
  dir.create(outdir)
}



groups=read.csv("12_gistic/12.gistic.sample.group.csv",stringsAsFactors = F)

##### broad level analysis:


broad_arm=read.table("12.sequenza.44.samples.gistic/broad_values_by_arm.txt",stringsAsFactors = F,sep="\t",header = T,check.names = F,row.names = 1)

pos=match(groups$Tumor,colnames(broad_arm))
groups=groups[pos,]

identical(groups$Tumor,colnames(broad_arm))


anno_col=groups[,c("group","ADT","Mutation_group")]
rownames(anno_col)=groups$Tumor


library(pheatmap)

pheatmap::pheatmap(broad_arm,annotation_col = anno_col,filename = file.path(outdir,"12.01.gistic.broad.arm.heatmap.pdf"),width = 12,height = 8)









### IDCP vs PCA

idcp.sample <- groups[groups$group == "IDCP", "Tumor"]
pca.sample <- groups[groups$group == "adenocarcinoma", "Tumor"]

idcp.broad=as.matrix(broad_arm[,idcp.sample])
pca.broad=as.matrix(broad_arm[,pca.sample])

cnv_fc=rowMeans(as.matrix(idcp.broad)) - rowMeans(as.matrix(pca.broad))


library(future.apply)
plan(multiprocess)
p_values <- future_lapply(seq(nrow(idcp.broad)), function(x){
  res <- wilcox.test(x = idcp.broad[x,], y =  pca.broad[x,],paired = T)
  res$p.value
})


p <- unlist(p_values)


diff_res <- data.frame(Symbol = rownames(idcp.broad),
                  Average_CNV_Number_Diff = cnv_fc,
                  pvalue = p,
                  logp=-log10(p))

write.csv(diff_res,file.path(outdir,"12.02.gistic.broad.arm.cnv.IDCP_vs_PCA.csv"),row.names = F,quote = F)

library(ggpubr)
library(ggrepel)

p1=ggscatter(diff_res,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (IDCP vs PCA)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(diff_res, logp >=-log10(0.05)), aes(label = Symbol))+ylim(0,3)
print(p1)
ggsave(file.path(outdir,"12.02.gistic.broad.arm.cnv.IDCP_vs_PCA.pdf"),width = 6,height = 6)




### ADT

ADT.sample <- groups[groups$ADT == "1", "Tumor"]
wild.sample <- groups[groups$ADT == "0", "Tumor"]

ADT.broad=as.matrix(broad_arm[,ADT.sample])
wild.broad=as.matrix(broad_arm[,wild.sample])

cnv_fc=rowMeans(as.matrix(ADT.broad)) - rowMeans(as.matrix(wild.broad))


library(future.apply)
plan(multiprocess)
p_values <- future_lapply(seq(nrow(ADT.broad)), function(x){
  res <- wilcox.test(x = ADT.broad[x,], y =  wild.broad[x,],paired = FALSE)
  res$p.value
})


p <- unlist(p_values)


ADT_diff_res <- data.frame(Symbol = rownames(ADT.broad),
                       Average_CNV_Number_Diff = cnv_fc,
                       pvalue = p,
                       logp=-log10(p))

write.csv(ADT_diff_res,file.path(outdir,"12.02.gistic.broad.arm.cnv.ADT_vs_noADT.csv"),row.names = F,quote = F)


library(ggpubr)
library(ggrepel)

p2=ggscatter(ADT_diff_res,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (ADT vs NO_ADT)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(ADT_diff_res, logp >=-log10(0.05)), aes(label = Symbol))+ylim(0,3)
print(p2)
ggsave(file.path(outdir,"12.02.gistic.broad.arm.cnv.ADT_vs_noADT.pdf"),width = 6,height = 6)


ggarrange(p1,p2,ncol=2,nrow=1,common.legend = T,labels = c("A)","B)"))
ggsave(file.path(outdir,"12.02.gistic.broad.arm.cnv.pdf"),width = 12,height = 6)






#### with IDCP

group_idcp=groups[groups$group=="IDCP",]

broad_arm.idcp=broad_arm[,group_idcp$Tumor]

## ADT
mut_sample <- group_idcp[group_idcp$ADT == "1", "Tumor"]
wild_sample <- group_idcp[group_idcp$ADT == "0", "Tumor"]

mut.broad=as.matrix(broad_arm.idcp[,mut_sample])
wild.broad=as.matrix(broad_arm.idcp[,wild_sample])

cnv_fc=rowMeans(as.matrix(mut.broad)) - rowMeans(as.matrix(wild.broad))


library(future.apply)
plan(multiprocess)
p_values <- future_lapply(seq(nrow(mut.broad)), function(x){
  res <- wilcox.test(x = mut.broad[x,], y =  wild.broad[x,],paired = FALSE)
  res$p.value
})


p <- unlist(p_values)


idcp.ADT_diff_res <- data.frame(Symbol = rownames(mut.broad),
                           Average_CNV_Number_Diff = cnv_fc,
                           pvalue = p,
                           logp=-log10(p))

write.csv(idcp.ADT_diff_res,file.path(outdir,"12.03.gistic.broad.arm.cnv.within_IDCP.ADT_vs_noADT.csv"),row.names = F,quote = F)



p_idcp_ADT=ggscatter(idcp.ADT_diff_res,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (ADT vs NO_ADT)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(idcp.ADT_diff_res, logp >=-log10(0.05)), aes(label = Symbol))+ylim(0,3)
print(p_idcp_ADT)
ggsave(file.path(outdir,"12.03.gistic.broad.arm.cnv.within_IDCP.ADT_vs_noADT.pdf"),width = 6,height = 6)





## A VS B
mut_sample <- group_idcp[group_idcp$Mutation_group == "A", "Tumor"]
wild_sample <- group_idcp[group_idcp$Mutation_group == "B", "Tumor"]

mut.broad=as.matrix(broad_arm.idcp[,mut_sample])
wild.broad=as.matrix(broad_arm.idcp[,wild_sample])

cnv_fc=rowMeans(as.matrix(mut.broad)) - rowMeans(as.matrix(wild.broad))


library(future.apply)
plan(multiprocess)
p_values <- future_lapply(seq(nrow(mut.broad)), function(x){
  res <- wilcox.test(x = mut.broad[x,], y =  wild.broad[x,],paired = FALSE)
  res$p.value
})


p <- unlist(p_values)


idcp.A_vs_B_diff_res <- data.frame(Symbol = rownames(mut.broad),
                                Average_CNV_Number_Diff = cnv_fc,
                                pvalue = p,
                                logp=-log10(p))

write.csv(idcp.A_vs_B_diff_res,file.path(outdir,"12.03.gistic.broad.arm.cnv.within_IDCP.A_vs_B.csv"),row.names = F,quote = F)



p_idcp_AB=ggscatter(idcp.A_vs_B_diff_res,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (A vs B)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(idcp.A_vs_B_diff_res, logp >=-log10(0.05)), aes(label = Symbol))+ylim(0,3)
print(p_idcp_AB)
ggsave(file.path(outdir,"12.03.gistic.broad.arm.cnv.within_IDCP.A_vs_B.pdf"),width = 6,height = 6)








#### within PCA

group_PCA=groups[groups$group=="adenocarcinoma",]

broad_arm.PCA=broad_arm[,group_PCA$Tumor]

## ADT
mut_sample <- group_PCA[group_PCA$ADT == "1", "Tumor"]
wild_sample <- group_PCA[group_PCA$ADT == "0", "Tumor"]

mut.broad=as.matrix(broad_arm.PCA[,mut_sample])
wild.broad=as.matrix(broad_arm.PCA[,wild_sample])

cnv_fc=rowMeans(as.matrix(mut.broad)) - rowMeans(as.matrix(wild.broad))


library(future.apply)
plan(multiprocess)
p_values <- future_lapply(seq(nrow(mut.broad)), function(x){
  res <- wilcox.test(x = mut.broad[x,], y =  wild.broad[x,],paired = FALSE)
  res$p.value
})


p <- unlist(p_values)


PCA.ADT_diff_res <- data.frame(Symbol = rownames(mut.broad),
                                Average_CNV_Number_Diff = cnv_fc,
                                pvalue = p,
                                logp=-log10(p))

write.csv(PCA.ADT_diff_res,file.path(outdir,"12.03.gistic.broad.arm.cnv.within_PCA.ADT_vs_noADT.csv"),row.names = F,quote = F)



p_pca_ADT=ggscatter(PCA.ADT_diff_res,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (ADT vs NO_ADT)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(PCA.ADT_diff_res, logp >=-log10(0.05)), aes(label = Symbol))+ylim(0,3)
print(p_pca_ADT)
ggsave(file.path(outdir,"12.03.gistic.broad.arm.cnv.within_PCA.ADT_vs_noADT.pdf"),width = 6,height = 6)





## A VS B
mut_sample <- group_PCA[group_PCA$Mutation_group == "A", "Tumor"]
wild_sample <- group_PCA[group_PCA$Mutation_group == "B", "Tumor"]

mut.broad=as.matrix(broad_arm.PCA[,mut_sample])
wild.broad=as.matrix(broad_arm.PCA[,wild_sample])

cnv_fc=rowMeans(as.matrix(mut.broad)) - rowMeans(as.matrix(wild.broad))


library(future.apply)
plan(multiprocess)
p_values <- future_lapply(seq(nrow(mut.broad)), function(x){
  res <- wilcox.test(x = mut.broad[x,], y =  wild.broad[x,],paired = FALSE)
  res$p.value
})


p <- unlist(p_values)


PCA.A_B_diff_res <- data.frame(Symbol = rownames(mut.broad),
                               Average_CNV_Number_Diff = cnv_fc,
                               pvalue = p,
                               logp=-log10(p))

write.csv(PCA.A_B_diff_res,file.path(outdir,"12.03.gistic.broad.arm.cnv.within_PCA.A_vs_B.csv"),row.names = F,quote = F)



p_pca_AB=ggscatter(PCA.A_B_diff_res,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (A vs B)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(PCA.A_B_diff_res, logp >=-log10(0.05)), aes(label = Symbol))+ylim(0,3)
print(p_pca_AB)
ggsave(file.path(outdir,"12.03.gistic.broad.arm.cnv.within_PCA.A_vs_B.pdf"),width = 6,height = 6)



ggarrange(p_idcp_ADT,p_idcp_AB,p_pca_ADT,p_pca_AB,ncol=2,nrow=2,common.legend = T,labels = c("A)","B)","C)","D)"))
ggsave(file.path(outdir,"12.03.gistic.broad.arm.cnv.group.compare.pdf"),width = 9,height = 9)








##### broad gene level analysis:


broad_arm_gene=read.table("12.sequenza.44.samples.gistic/broad_data_by_genes.txt",stringsAsFactors = F,sep="\t",header = T,check.names = F,row.names = 1)


broad_arm_gene_anno=broad_arm_gene[,c(1,2)]
broad_arm_gene_anno$gene=rownames(broad_arm_gene_anno)


write.csv(broad_arm_gene_anno,file.path(outdir,"12.04.gistic.broad.arm.gene.list.csv"),row.names = F,quote = F)


broad_arm_gene$`Gene ID`=NULL
broad_arm_gene$Cytoband=NULL



pos=match(colnames(broad_arm_gene),groups$Tumor)
groups=groups[pos,]

identical(groups$Tumor,colnames(broad_arm_gene))


anno_col=groups[,c("group","ADT","Mutation_group")]
rownames(anno_col)=groups$Tumor


## ADT
mut_sample <- groups[groups$ADT == "1", "Tumor"]
wild_sample <- groups[groups$ADT == "0", "Tumor"]

mut.broad=as.matrix(broad_arm_gene[,mut_sample])
wild.broad=as.matrix(broad_arm_gene[,wild_sample])

cnv_fc=rowMeans(as.matrix(mut.broad)) - rowMeans(as.matrix(wild.broad))

library(future.apply)
plan(multiprocess)
p_values <- future_lapply(seq(nrow(mut.broad)), function(x){
  res <- wilcox.test(x = mut.broad[x,], y =  wild.broad[x,],paired = FALSE)
  res$p.value
})


p <- unlist(p_values)

ADT_diff_res <- data.frame(Symbol = rownames(mut.broad),
                               Average_CNV_Number_Diff = cnv_fc,
                               pvalue = p,
                               logp=-log10(p))



## IDCP VS pca
mut_sample <- groups[groups$group == "IDCP", "Tumor"]
wild_sample <- groups[groups$group == "adenocarcinoma", "Tumor"]

mut.broad=as.matrix(broad_arm_gene[,mut_sample])
wild.broad=as.matrix(broad_arm_gene[,wild_sample])

cnv_fc=rowMeans(as.matrix(mut.broad)) - rowMeans(as.matrix(wild.broad))

library(future.apply)
plan(multiprocess)
p_values <- future_lapply(seq(nrow(mut.broad)), function(x){
  res <- wilcox.test(x = mut.broad[x,], y =  wild.broad[x,],paired = FALSE)
  res$p.value
})


p <- unlist(p_values)

IDCP_PCA_diff_res <- data.frame(Symbol = rownames(mut.broad),
                           Average_CNV_Number_Diff = cnv_fc,
                           pvalue = p,
                           logp=-log10(p))

table(ADT_diff_res$pvalue<0.05)

table(IDCP_PCA_diff_res$pvalue<0.05)

gene_selected=union(ADT_diff_res$Symbol[ADT_diff_res$pvalue<0.05],IDCP_PCA_diff_res$Symbol[IDCP_PCA_diff_res$pvalue<0.05])


broad_arm_gene_plot=broad_arm_gene[gene_selected,]

anno_row=broad_arm_gene_anno[gene_selected,]
anno_row$`Gene ID`=NULL
anno_row$gene=NULL

anno_col$ADT=as.character(anno_col$ADT)
anno_col$Mutation_group

pheatmap::pheatmap(broad_arm_gene_plot,annotation_col = anno_col,cluster_rows = FALSE,
                   filename = file.path(outdir,"12.04.gistic.broad.arm.gene.heatmap.pdf"),width = 12,height = 8,show_rownames = F)









