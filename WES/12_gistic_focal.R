#  setwd("../")


rm(list=ls())

outdir="12_gistic_focal"

if(!dir.exists(outdir)){
  dir.create(outdir)
}



groups=read.csv("12_gistic/12.gistic.sample.group.csv",stringsAsFactors = F)

##### focal gene level analysis:


focal_genes=read.table("12.sequenza.44.samples.gistic/all_thresholded.by_genes.txt",stringsAsFactors = F,sep="\t",header = T,check.names = F,row.names = 1)

focal_genes_anno=focal_genes[,1:2]
focal_genes_anno$gene_symbol=rownames(focal_genes)



samples_used=intersect(groups$Tumor,colnames(focal_genes))

### final data
focal_genes=focal_genes[,samples_used]
groups=groups[match(samples_used,groups$Tumor),]

identical(colnames(focal_genes),groups$Tumor)



### IDCP vs PCA

## 
mut_sample <- groups[groups$group == "IDCP", "Tumor"]
wild_sample <- groups[groups$group == "adenocarcinoma", "Tumor"]



#### function 
cnvdiff=function(focal_genes,mut_sample,wild_sample,isPaired){
  
  mut=as.matrix(focal_genes[,mut_sample])
  wild=as.matrix(focal_genes[,wild_sample])
  
  cnv_fc=rowMeans(as.matrix(mut)) - rowMeans(as.matrix(wild))
  
  library(future.apply)
  plan(multiprocess)
  p_values <- future_lapply(seq(nrow(mut)), function(x){
    res <- wilcox.test(x = mut[x,], y =  wild[x,],paired = isPaired)
    res$p.value
  })
  
  
  p <- unlist(p_values)
  padj=p.adjust(p,method = "fdr")
  
  
  df <- data.frame(Symbol = rownames(mut),
                                  Average_CNV_Number_Diff = cnv_fc,
                                  pvalue = p,
                                  padj=padj,
                                  logp=-log10(p),stringsAsFactors = F)
  
  return(df)
}
#####


idcp_pca=cnvdiff(focal_genes,mut_sample,wild_sample,TRUE)
idcp_pca$significance="No"
idcp_pca$significance[idcp_pca$pvalue<0.05]="Yes"

idcp_pca$padj_significance="No"
idcp_pca$padj_significance[idcp_pca$padj<0.05]="Yes"


write.csv(idcp_pca,file.path(outdir,"12.01.gistic.gene.cnv.idcp_vs_pca.csv"),row.names = F,quote = F)



library(pheatmap)


idcp_pca_gene_df=focal_genes[idcp_pca$Symbol[idcp_pca$significance=="Yes"],]

anno_col=data.frame(Group=groups$group,
                    ADT=as.character(groups$ADT),
                    Mutation_Group=groups$Mutation_group,stringsAsFactors = F)
rownames(anno_col)=groups$Tumor



pheatmap::pheatmap(idcp_pca_gene_df,annotation_col = anno_col,show_rownames = F,
                   filename = file.path(outdir,"12.01.gistic.gene.cnv.idcp_vs_pca.pdf"),width = 9,height = 9)





##### KEGG and GO enrichment

library(dplyr)
library(magrittr)
library(maftools)
library(edgeR)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(openxlsx)
library(samr)
library(ggplot2)


# 背景基因
dea.df <- bitr(unique(rownames(focal_genes)), fromType = "SYMBOL", ###输入只能说entrez id，所以需要iD转换
               toType = c("SYMBOL","ENTREZID"),
               OrgDb = org.Hs.eg.db) 


deg.df <- bitr(unique(rownames(idcp_pca_gene_df)), fromType = "SYMBOL", ###输入只能说entrez id，所以需要iD转换
               toType = c("SYMBOL","ENTREZID"),
               OrgDb = org.Hs.eg.db) 



deg.go <- enrichGO(deg.df$ENTREZID, 
                   universe = dea.df$ENTREZID,
                   OrgDb = "org.Hs.eg.db", 
                   ont="all",
                   readable = T,
                   qvalueCutoff = 0.05) 



pdf(file = file.path(outdir,"12.01.gistic.gene.cnv.idcp_vs_pca.GO.pdf"),onefile = F,width = 9,height = 9)
enrichplot::dotplot(deg.go, showCategory = 30) + 
  scale_y_discrete(labels=function(y) stringr::str_wrap(y, width=40))
dev.off()



deg.go.df=as.data.frame(deg.go)
write.csv(deg.go.df,file.path(outdir,"12.01.gistic.gene.cnv.idcp_vs_pca.GO.csv"),quote = F)



kk <- enrichKEGG(gene         = deg.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)



pdf(file = file.path(outdir,"12.01.gistic.gene.cnv.idcp_vs_pca.KEGG.pdf"),onefile = F,width = 8,height = 9)
enrichplot::dotplot(kk, showCategory = 30) + 
  scale_y_discrete(labels=function(y) stringr::str_wrap(y, width=40))
dev.off()

kkread <- setReadable(kk, org.Hs.eg.db, keyType = "ENTREZID")

kk.df=as.data.frame(kkread)
write.csv(kk.df,file.path(outdir,"12.01.gistic.gene.cnv.idcp_vs_pca.KEGG.csv"),quote = F)





### ADT vs NONADT

## 
mut_sample <- groups[groups$ADT == "1", "Tumor"]
wild_sample <- groups[groups$ADT == "0", "Tumor"]

gene_ADT=cnvdiff(focal_genes,mut_sample,wild_sample,FALSE)

gene_ADT$significance="No"
gene_ADT$significance[gene_ADT$pvalue<0.05]="Yes"

gene_ADT$padj_significance="No"
gene_ADT$padj_significance[gene_ADT$padj<0.05]="Yes"

write.csv(gene_ADT,file.path(outdir,"12.02.gistic.gene.cnv.ADT_vs_nonADT.csv"),row.names = F,quote = F)



gene_ADT_deg_plot=focal_genes[gene_ADT$Symbol[gene_ADT$significance=="Yes"],]

pheatmap::pheatmap(gene_ADT_deg_plot,annotation_col = anno_col,show_rownames = FALSE,
                   filename = file.path(outdir,"12.02.gistic.gene.cnv.ADT_vs_nonADT.pdf"),width = 8,height = 8)



#### enrichment

dea.df <- bitr(unique(rownames(focal_genes)), fromType = "SYMBOL", ###输入只能说entrez id，所以需要iD转换
               toType = c("SYMBOL","ENTREZID"),
               OrgDb = org.Hs.eg.db) 


deg.df <- bitr(unique(rownames(gene_ADT_deg_plot)), fromType = "SYMBOL", ###输入只能说entrez id，所以需要iD转换
               toType = c("SYMBOL","ENTREZID"),
               OrgDb = org.Hs.eg.db) 



deg.go <- enrichGO(deg.df$ENTREZID, 
                   universe = dea.df$ENTREZID,
                   OrgDb = "org.Hs.eg.db", 
                   ont="all",
                   readable = T,
                   qvalueCutoff = 0.05) 



pdf(file = file.path(outdir,"12.02.gistic.gene.cnv.ADT_vs_nonADT.GO.pdf"),onefile = F,width = 9,height = 9)
enrichplot::dotplot(deg.go, showCategory = 30) + 
  scale_y_discrete(labels=function(y) stringr::str_wrap(y, width=40))
dev.off()



deg.go.df=as.data.frame(deg.go)
write.csv(deg.go.df,file.path(outdir,"12.02.gistic.gene.cnv.ADT_vs_nonADT.GO.csv"),quote = F)



kk <- enrichKEGG(gene         = deg.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)



pdf(file = file.path(outdir,"12.02.gistic.gene.cnv.ADT_vs_nonADT.KEGG.pdf"),onefile = F,width = 8,height = 9)
enrichplot::dotplot(kk, showCategory = 30) + 
  scale_y_discrete(labels=function(y) stringr::str_wrap(y, width=40))
dev.off()

kkread <- setReadable(kk, org.Hs.eg.db, keyType = "ENTREZID")

kk.df=as.data.frame(kkread)
write.csv(kk.df,file.path(outdir,"12.02.gistic.gene.cnv.ADT_vs_nonADT.KEGG.csv"),quote = F)







venn.diagram(x=list(IDCP_vs_PCA=rownames(idcp_pca_gene_df),
                    ADT_vs_NonADT=rownames(gene_ADT_deg_plot)),
             filename = file.path(outdir,"12.02.gistic.gene.cnv.IDCP_PCA.ADT_nonADT.KEGG.tiff"),
             
             fill=c("red","yellow"),margin=c(0.1,0.1,0.1,0.1)
)










### within IDCP

idcp=focal_genes[,groups$Tumor[groups$group=="IDCP"]]

idcp_group=groups[groups$group=="IDCP",]


##  ADT VS non_ADT
mut_sample <- idcp_group[idcp_group$ADT == "1", "Tumor"]
wild_sample <- idcp_group[idcp_group$ADT == "0", "Tumor"]

focal_gene_IDCP_ADT=cnvdiff(idcp,mut_sample,wild_sample,FALSE)

focal_gene_IDCP_ADT$significance="No"
focal_gene_IDCP_ADT$significance[focal_gene_IDCP_ADT$pvalue<0.05]="Yes"

focal_gene_IDCP_ADT$padj_significance="No"
focal_gene_IDCP_ADT$padj_significance[focal_gene_IDCP_ADT$padj<0.05]="Yes"


write.csv(focal_gene_IDCP_ADT,file.path(outdir,"12.03.gistic.focal.gene.cnv.withinIDCP.ADT_vs_nonADT.csv"),row.names = F,quote = F)


focal_gene_IDCP_ADT_plot=focal_genes[focal_gene_IDCP_ADT$Symbol[focal_gene_IDCP_ADT$significance=="Yes"],]

pheatmap::pheatmap(focal_gene_IDCP_ADT_plot,annotation_col = anno_col,cluster_rows = FALSE,show_rownames = FALSE,
                   filename = file.path(outdir,"12.03.gistic.focal.gene.cnv.withinIDCP.ADT_vs_nonADT.pdf"),width = 8,height = 8)



###
##  A VS B
mut_sample <- idcp_group[idcp_group$Mutation_group == "A", "Tumor"]
wild_sample <- idcp_group[idcp_group$Mutation_group == "B", "Tumor"]

focal_gene_IDCP_A_B=cnvdiff(idcp,mut_sample,wild_sample,FALSE)

focal_gene_IDCP_A_B$significance="No"
focal_gene_IDCP_A_B$significance[focal_gene_IDCP_A_B$pvalue<0.05]="Yes"

focal_gene_IDCP_A_B$padj_significance="No"
focal_gene_IDCP_A_B$padj_significance[focal_gene_IDCP_A_B$padj<0.05]="Yes"

write.csv(focal_gene_IDCP_A_B,file.path(outdir,"12.03.gistic.focal.gene.cnv.withinIDCP.A_vs_B.csv"),row.names = F,quote = F)


focal_gene_IDCP_A_B_plot=focal_genes[focal_gene_IDCP_A_B$Symbol[focal_gene_IDCP_A_B$significance=="Yes"],]


pheatmap::pheatmap(focal_gene_IDCP_A_B_plot,annotation_col = anno_col,cluster_rows = FALSE,show_rownames = FALSE,
                   filename = file.path(outdir,"12.03.gistic.focal.gene.cnv.withinIDCP.A_vs_B.pdf"),width = 8,height = 8)





### within PCA

pca=focal_genes[,groups$Tumor[groups$group=="adenocarcinoma"]]

pca_group=groups[groups$group=="adenocarcinoma",]


##  ADT VS non_ADT
mut_sample <- pca_group[pca_group$ADT == "1", "Tumor"]
wild_sample <- pca_group[pca_group$ADT == "0", "Tumor"]

focal_gene_pca_ADT=cnvdiff(pca,mut_sample,wild_sample,FALSE)

focal_gene_pca_ADT$significance="No"
focal_gene_pca_ADT$significance[focal_gene_pca_ADT$pvalue<0.05]="Yes"

focal_gene_pca_ADT$padj_significance="No"
focal_gene_pca_ADT$padj_significance[focal_gene_pca_ADT$padj<0.05]="Yes"



write.csv(focal_gene_pca_ADT,file.path(outdir,"12.04.gistic.focal.gene.cnv.withinPCA.ADT_vs_nonADT.csv"),row.names = F,quote = F)


focal_gene_pca_ADT_plot=focal_genes[focal_gene_pca_ADT$Symbol[focal_gene_pca_ADT$significance=="Yes"],]


pheatmap::pheatmap(focal_gene_pca_ADT_plot,annotation_col = anno_col,cluster_rows = FALSE,show_rownames = FALSE,
                   filename = file.path(outdir,"12.04.gistic.focal.gene.cnv.withinPCA.ADT_vs_nonADT.pdf"),width = 8,height = 8)



##  A VS B
mut_sample <- pca_group[pca_group$Mutation_group == "A", "Tumor"]
wild_sample <- pca_group[pca_group$Mutation_group == "B", "Tumor"]

focal_gene_pca_AB=cnvdiff(pca,mut_sample,wild_sample,FALSE)

focal_gene_pca_AB$significance="No"
focal_gene_pca_AB$significance[focal_gene_pca_AB$pvalue<0.05]="Yes"

focal_gene_pca_AB$padj_significance="No"
focal_gene_pca_AB$padj_significance[focal_gene_pca_AB$padj<0.05]="Yes"


write.csv(focal_gene_pca_AB,file.path(outdir,"12.04.gistic.focal.gene.cnv.withinPCA.A_vs_B.csv"),row.names = F,quote = F)


focal_gene_pca_AB_plot=focal_genes[focal_gene_pca_AB$Symbol[focal_gene_pca_AB$significance=="Yes"],]


pheatmap::pheatmap(focal_gene_pca_AB_plot,annotation_col = anno_col,cluster_rows = FALSE,show_rownames = FALSE,
                   filename = file.path(outdir,"12.04.gistic.focal.gene.cnv.withinPCA.A_vs_B.pdf"),width = 8,height = 8)


library(VennDiagram)
library(RColorBrewer)

venn.diagram(x=list(IDCP_ADT=row.names(focal_gene_IDCP_ADT_plot),IDCP_AB=rownames(focal_gene_IDCP_A_B_plot),
                    PCA_ADT=rownames(focal_gene_pca_ADT_plot),PCA_AB=rownames(focal_gene_pca_AB_plot)),
             filename = file.path(outdir,"12.05.gistic.focal.gene.cnv.venn.tiff"),
             fill = brewer.pal(4, "Set2"),
             col = brewer.pal(4, "Set3"),
             )



venn.diagram(x=list(IDCP_AB=rownames(focal_gene_IDCP_A_B_plot),
                    PCA_AB=rownames(focal_gene_pca_AB_plot)),
             filename = file.path(outdir,"12.05.gistic.focal.gene.cnv.venn.AB.tiff"),
             fill=c("red","yellow"),margin=c(0.1,0.1,0.1,0.1)
)


venn.diagram(x=list(IDCP_ADT=row.names(focal_gene_IDCP_ADT_plot),PCA_ADT=rownames(focal_gene_pca_ADT_plot)),
             filename = file.path(outdir,"12.05.gistic.focal.gene.cnv.venn.ADT.tiff"),
             fill=c("red","yellow"),margin=c(0.1,0.1,0.1,0.1)
)








############## fisher exact test


mut_sample <- groups[groups$group == "IDCP", "Tumor"]
wild_sample <- groups[groups$group == "adenocarcinoma", "Tumor"]




#### function 
fisherresult=function(focal_genes,mut_sample,wild_sample){
  
  mut=as.matrix(focal_genes[,mut_sample])
  wild=as.matrix(focal_genes[,wild_sample])
  
  cnv_fc=rowMeans(as.matrix(mut)) - rowMeans(as.matrix(wild))
  
  library(future.apply)
  plan(multiprocess)
  p_values <- future_lapply(seq(nrow(mut)), function(x){
    
    #res <- wilcox.test(x = mut[x,], y =  wild[x,],paired = isPaired)
    mut_cnv=mut[x,]
    mut_gain=length(which(mut_cnv>0))
    mut_loss=length(which(mut_cnv<0))
    
    wild_cnv=wild[x,]
    
    wild_gain=length(which(wild_cnv>0))
    wild_loss=length(which(wild_cnv<0))
    
    mat=cbind(c(mut_gain,mut_loss),c(wild_gain,wild_loss))
    colnames(mat)=c("Mut","Wild")
    rownames(mat)=c("Gain","Loss")
    
    res=fisher.test(mat)
    return(res$p.value)
    
  })
  
  
  p <- unlist(p_values)
  padj=p.adjust(p,method = "fdr")
  
  
  df <- data.frame(Symbol = rownames(mut),
                   Average_CNV_Number_Diff = cnv_fc,
                   pvalue = p,
                   padj=padj,
                   logp=-log10(p),stringsAsFactors = F)
  
  
  df$significance="No"
  df$significance[df$pvalue<0.05]="Yes"
  
  df$padj_significance="No"
  df$padj_significance[df$padj<0.05]="Yes"
  
  
  
  return(df)
}
#####

IDCP_PCF_fisher=fisherresult(focal_genes,mut_sample,wild_sample)



write.csv(IDCP_PCF_fisher,file.path(outdir,"12.06.gistic.focal.gene.cnv.IDCP_vs_PCA.fisher.csv"),row.names = F,quote = F)




venn.diagram(x=list(Wilcox=idcp_pca$Symbol[idcp_pca$significance=="Yes"],
                    Fisher=IDCP_PCF_fisher$Symbol[IDCP_PCF_fisher$significance=="Yes"]),
             filename = file.path(outdir,"12.07.gistic.focal.gene.cnv.wilcox_vs_fisher.venn.tiff"),
             
             fill=c("red","yellow"),margin=c(0.1,0.1,0.1,0.1)
)




