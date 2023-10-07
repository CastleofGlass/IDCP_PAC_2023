#  setwd("../")


rm(list=ls())

outdir="12_gistic_peak"

if(!dir.exists(outdir)){
  dir.create(outdir)
}



groups=read.csv("12_gistic/12.gistic.sample.group.csv",stringsAsFactors = F)

##### focal gene level analysis:


peak=read.table("12.sequenza.44.samples.gistic/all_lesions.conf_90.txt",stringsAsFactors = F,sep="\t",header = T,check.names = F)

pos=grep("Actual Copy Change Given",peak$`Amplitude Threshold`)


peak=peak[pos,]


peak_anno=peak[,1:9]
peak_anno$group=sapply(peak_anno$`Unique Name`,function(x)unlist(strsplit(x," Peak"))[1])


peak_cnv=peak[,10:ncol(peak)]
rownames(peak_cnv)=peak_anno$`Unique Name`


samples_used=intersect(groups$Tumor,colnames(peak_cnv))

peak_cnv=peak_cnv[,samples_used]
groups=groups[match(samples_used,groups$Tumor),]


peak_col_anno=data.frame(Group=groups$group)
peak_col_anno$ADT=as.character(groups$ADT)
peak_col_anno$Mutation_group=as.character(groups$Mutation_group)
rownames(peak_col_anno)=groups$Tumor


peak_anno_row=peak_anno[match(rownames(peak_cnv),peak_anno$`Unique Name`),]
peak_row_anno=peak_anno_row[,c("q values","group")]
rownames(peak_row_anno)=peak_anno_row$`Unique Name`

library(pheatmap)

pheatmap::pheatmap(as.matrix(peak_cnv),annotation_row = peak_row_anno, annotation_col = peak_col_anno)#,
                   filename = file.path(outdir,"12.01.gistic.peak.heatmap.pdf"),width = 9,height = 8,show_rownames = F)





### IDCP vs PCA

## 
mut_sample <- groups[groups$group == "IDCP", "Tumor"]
wild_sample <- groups[groups$group == "adenocarcinoma", "Tumor"]



#### function 
cnvdiff=function(dat,mut_sample,wild_sample,isPaired){
  
  mut=as.matrix(dat[,mut_sample])
  wild=as.matrix(dat[,wild_sample])
  
  cnv_fc=rowMeans(as.matrix(mut)) - rowMeans(as.matrix(wild))
  
  library(future.apply)
  plan(multiprocess)
  p_values <- future_lapply(seq(nrow(mut)), function(x){
    res <- wilcox.test(x = mut[x,], y =  wild[x,],paired = isPaired)
    res$p.value
  })
  
  
  p <- unlist(p_values)
  
  
  df <- data.frame(Symbol = rownames(mut),
                                  Average_CNV_Number_Diff = cnv_fc,
                                  pvalue = p,
                                  logp=-log10(p),stringsAsFactors = F)
  
  return(df)
}
#####


peak_idcp_pca=cnvdiff(peak_cnv,mut_sample,wild_sample,TRUE)

write.csv(peak_idcp_pca,file.path(outdir,"12.02.gistic.peak.IDCP_vs_PCA.csv"),row.names = F,quote = F)


library(ggpubr)
library(ggrepel)

peak_idcp_pca$Symbol=gsub(" - CN values","",peak_idcp_pca$Symbol)

p1=ggscatter(peak_idcp_pca,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (IDCP vs PCA)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(peak_idcp_pca, logp >=-log10(0.05)), aes(label = Symbol))
print(p1)
ggsave(file.path(outdir,"12.02.gistic.peak.IDCP_vs_PCA.pdf"),width = 6,height = 6)




### ADT vs nonADT
## 
mut_sample <- groups[groups$ADT == "1", "Tumor"]
wild_sample <- groups[groups$ADT == "0", "Tumor"]

peak_ADT_nonADT=cnvdiff(peak_cnv,mut_sample,wild_sample,FALSE)

write.csv(peak_ADT_nonADT,file.path(outdir,"12.03.gistic.peak.ADT_vs_nonADT.csv"),row.names = F,quote = F)


peak_ADT_nonADT$Symbol=gsub(" - CN values","",peak_ADT_nonADT$Symbol)

p2=ggscatter(peak_ADT_nonADT,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (ADT vs NonADT)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(peak_ADT_nonADT, logp >=-log10(0.05)), aes(label = Symbol))
print(p2)
ggsave(file.path(outdir,"12.03.gistic.peak.ADT_vs_nonADT.pdf"),width = 6,height = 6)



ggarrange(p1+ylim(0,3.2),p2+ylim(0,3.2),ncol=2,nrow=1,common.legend = T,labels = c("A)","B)"))
ggsave(file.path(outdir,"12.03.gistic.peak.all.cnv.pdf"),width = 10,height = 9)



### withinIDCP
## 

idcp_group=groups[groups$group=="IDCP",]

idcp=peak_cnv[,idcp_group$Tumor]

##ADT vs nonADT

mut_sample <- idcp_group[idcp_group$ADT == "1", "Tumor"]
wild_sample <- idcp_group[idcp_group$ADT == "0", "Tumor"]

idcp_ADT_nonADT=cnvdiff(idcp,mut_sample,wild_sample,FALSE)

write.csv(idcp_ADT_nonADT,file.path(outdir,"12.04.gistic.IDCP.ADT_vs_nonADT.csv"),row.names = F,quote = F)

idcp_ADT_nonADT$Symbol=gsub(" - CN values","",idcp_ADT_nonADT$Symbol)


p3=ggscatter(idcp_ADT_nonADT,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (ADT vs NonADT)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(idcp_ADT_nonADT, logp >=-log10(0.05)), aes(label = Symbol))
print(p3)
ggsave(file.path(outdir,"12.04.gistic.IDCP.ADT_vs_nonADT.pdf"),width = 6,height = 6)



## A vs B

mut_sample <- idcp_group[idcp_group$Mutation_group == "A", "Tumor"]
wild_sample <- idcp_group[idcp_group$Mutation_group == "B", "Tumor"]

idcp_AB=cnvdiff(idcp,mut_sample,wild_sample,FALSE)

write.csv(idcp_AB,file.path(outdir,"12.05.gistic.IDCP.A_vs_B.csv"),row.names = F,quote = F)


idcp_AB$Symbol=gsub(" - CN values","",idcp_AB$Symbol)

p4=ggscatter(idcp_AB,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (A vs B)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(idcp_AB, logp >=-log10(0.05)), aes(label = Symbol))
print(p4)
ggsave(file.path(outdir,"12.05.gistic.IDCP.A_vs_B.pdf"),width = 6,height = 6)





### within PCA

## ADT vs nonADT

pca_group=groups[groups$group=="adenocarcinoma",]

pca=peak_cnv[,pca_group$Tumor]

##ADT vs nonADT

mut_sample <- pca_group[pca_group$ADT == "1", "Tumor"]
wild_sample <- pca_group[pca_group$ADT == "0", "Tumor"]

pca_ADT_nonADT=cnvdiff(pca,mut_sample,wild_sample,FALSE)

write.csv(pca_ADT_nonADT,file.path(outdir,"12.06.gistic.PCA.ADT_vs_nonADT.csv"),row.names = F,quote = F)


pca_ADT_nonADT$Symbol=gsub(" - CN values","",pca_ADT_nonADT$Symbol)

p5=ggscatter(pca_ADT_nonADT,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (ADT vs NonADT)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(pca_ADT_nonADT, logp >=-log10(0.05)), aes(label = Symbol))
print(p5)
ggsave(file.path(outdir,"12.06.gistic.PCA.ADT_vs_nonADT.pdf"),width = 6,height = 6)




##A vs B

mut_sample <- pca_group[pca_group$Mutation_group == "A", "Tumor"]
wild_sample <- pca_group[pca_group$Mutation_group == "B", "Tumor"]

pca_AB_nonADT=cnvdiff(pca,mut_sample,wild_sample,FALSE)

write.csv(pca_AB_nonADT,file.path(outdir,"12.07.gistic.PCA.A_vs_B.csv"),row.names = F,quote = F)


pca_AB_nonADT$Symbol=gsub(" - CN values","",pca_AB_nonADT$Symbol)

p6=ggscatter(pca_AB_nonADT,x="Average_CNV_Number_Diff",y="logp",shape = 21,size = 3)+ylab("-Log10 P value (wilcox)")+xlab("Average_CNV_Number_Diff (A vs B)")+
  geom_hline(yintercept = -log10(0.05),lty=4)+
  geom_text_repel(data = subset(pca_AB_nonADT, logp >=-log10(0.05)), aes(label = Symbol))
print(p6)
ggsave(file.path(outdir,"12.07.gistic.PCA.A_vs_B.pdf"),width = 6,height = 6)


ggarrange(p3+ylim(0,2.6),p4+ylim(0,2.6),p5+ylim(0,2.6),p6+ylim(0,2.6),ncol=2,nrow=2,common.legend = T,labels = c("A)","B)","C)","D)"))
ggsave(file.path(outdir,"12.08.gistic.peak.idcp.pca.cnv.pdf"),width = 12,height = 12)





