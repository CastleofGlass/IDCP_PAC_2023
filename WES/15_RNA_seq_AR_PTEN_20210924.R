rm(list=ls())
options(stringsAsFactors = F)

##  setwd("../")
result_dir="15_RNA_AR_PTEN"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}


RNA=read.csv("data/gene_tpm_anno_HXID.csv",stringsAsFactors = F,check.names = F)

RNA$gene_id=NULL


library(reshape)

rna_df=melt(RNA,id.vars = "entreID")

colnames(rna_df)=c("entreID","Sample","TPM")
rna_df$id=sapply(as.character(rna_df$Sample),function(x)paste(unlist(strsplit(x,"-"))[1:2],collapse = "-"))
rna_df$group=sapply(as.character(rna_df$Sample),function(x)unlist(strsplit(x,"-"))[3])
rna_df$log2RPM=log2(rna_df$TPM+1)


library(ggpubr)


p<-ggboxplot(rna_df,x="Sample",y="TPM",color="group",palette="jco")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
print(p)
ggsave(file.path(result_dir,"01.RNA.TPM.boxplot.pdf"),width = 12,height = 10)



p<-ggboxplot(rna_df,x="Sample",y="log2RPM",color="group",palette="jco")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
print(p)
ggsave(file.path(result_dir,"02.RNA.log2TPM.boxplot.pdf"),width = 12,height = 10)



AR=rna_df[which(rna_df$entreID==367),]
PTEN=rna_df[which(rna_df$entreID==5728),]

p<-ggline(AR,x="group",y="TPM",title="AR")+facet_wrap("id")
print(p)
ggsave(file.path(result_dir,"03.RNA.AR.expression.pdf"),width = 12,height = 10)


p<-ggline(AR,x="group",y="log2RPM",title="AR")+facet_wrap("id")
print(p)
ggsave(file.path(result_dir,"03.RNA.log2AR.expression.pdf"),width = 12,height = 10)



p<-ggline(PTEN,x="group",y="TPM",title="PTEN")+facet_wrap("id")
print(p)
ggsave(file.path(result_dir,"04.RNA.PTEN.expression.pdf"),width = 12,height = 10)


p<-ggline(PTEN,x="group",y="log2RPM",title="PTEN")+facet_wrap("id")
print(p)
ggsave(file.path(result_dir,"04.RNA.log2PTEN.expression.pdf"),width = 12,height = 10)





