rm(list=ls())
options(stringsAsFactors = F)

#setwd("../")

args=c("07_sequenza/cnv_anno","data/idcp.pairs.with.HXid.csv")

samples=read.csv(args[2],stringsAsFactors = F)
samples

cnv_files=list.files(path=args[1],pattern = "filtered.tsv$",all.files = T,full.names = T,recursive = T)
cnv_files

cnv_files_df=data.frame(file=cnv_files,stringsAsFactors = F)
cnv_files_df

cnv_files_df$prefix=sapply(cnv_files_df$file, function(x)gsub(".sequenza.result.anno.bed.filtered.tsv","",basename(x)))
cnv_files_df

result=c()

for(i in 1:nrow(samples)){
  print(i)
  idcp_id=samples$idcp[i]
  ade_id=samples$adenocarcinoma[i]
  
  idcp_file=cnv_files_df$file[cnv_files_df$prefix==idcp_id]
  ade_file=cnv_files_df$file[cnv_files_df$prefix==ade_id]
  
  idcp_df=read.table(idcp_file,stringsAsFactors = F,header = T,sep="\t")
  idcp_df$full_name=paste(idcp_df$Gene,idcp_df$CNV_classification,sep="_")
  idcp_df=idcp_df[idcp_df$CNV_classification!=".",]
  
  ade_df=read.table(ade_file,stringsAsFactors = F,header = T,sep="\t")
  ade_df$full_name=paste(ade_df$Gene,ade_df$CNV_classification,sep="_")
  ade_df=ade_df[ade_df$CNV_classification!=".",]
  
  
  idcp_gain=idcp_df[idcp_df$CNV_classification=="gain",]
  idcp_loss=idcp_df[idcp_df$CNV_classification=="loss",]
  
  ade_gain=ade_df[ade_df$CNV_classification=="gain",]
  ade_loss=ade_df[ade_df$CNV_classification=="loss",]
  
  
  res=c(idcp_id,nrow(idcp_df),nrow(idcp_gain),nrow(idcp_loss),ade_id,nrow(ade_df),nrow(ade_gain),nrow(ade_loss),
        length(union(idcp_df$full_name,ade_df$full_name)),length(intersect(idcp_df$full_name,ade_df$full_name)),
        length(union(idcp_gain$full_name,ade_gain$full_name)),length(intersect(idcp_gain$full_name,ade_gain$full_name)),
        length(union(idcp_loss$full_name,ade_loss$full_name)),length(intersect(idcp_loss$full_name,ade_loss$full_name))  )
  
  result=rbind(result,res)
}


result=as.data.frame(result)
colnames(result)=c("IDCP_ID","IDCP_Number_of_totalCNV_gene_level","IDCP_Number_of_GAIN_gene_level","IDCP_Number_of_LOSS_gene_level",
                   "PCA_ID","PCA_Number_of_totalCNV_gene_level","PCA_Number_of_GAIN_gene_level","PCA_Number_of_LOSS_gene_level",
                   "Total_CNV_per_sample","Total_CNV_Overlapped",
                   "Gain_CNV_per_sample","Gain_CNV_Overlapped",
                   "LOSS_CNV_per_sample","Loss_CNV_overlapped")



outdir="07_sequenza"

write.csv(result,file.path(outdir,"sequenza_gene_level_cnv_overlapped.csv"),quote = F,row.names = F)



df=read.csv(file.path(outdir,"sequenza_gene_level_cnv_overlapped.csv"),stringsAsFactors = F,header = T)

df$TOtal_CNV_percentage=100*df$Total_CNV_Overlapped/df$Total_CNV_per_sample
df$GAIN_CNV_percentage=100*df$Gain_CNV_Overlapped/df$Gain_CNV_per_sample
df$LOSS_CNV_percentage=100*df$Loss_CNV_overlapped/df$LOSS_CNV_per_sample


pos=match(df$IDCP_ID,samples$idcp)
df=cbind(df,samples[pos,])


write.csv(df,file.path(outdir,"sequenza_gene_level_cnv_overlapped_update.csv"),quote = F,row.names = F)



#### total

df_plot=df[,c(2,6,22)]
colnames(df_plot)
colnames(df_plot)[1:2]=c("IDCP","PCA")


library(reshape)
df_plot_dat=melt(df_plot,id.vars = "hxid")


my_comparisons <- list( c("IDCP", "PCA"))

library(ggpubr)

ggpaired(df_plot_dat, x = "variable", y = "value",shape = "variable",color = "variable",id="hxid",palette = "npg",label="hxid",
         label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Total Gene CNV numbers")

ggsave(file.path(outdir,paste("sequenza","total","gene_level_cnv.pdf",sep=".")),width = 6,height = 6)





#### gain

df_plot=df[,c(3,7,22)]
colnames(df_plot)
colnames(df_plot)[1:2]=c("IDCP","PCA")


library(reshape)
df_plot_dat=melt(df_plot,id.vars = "hxid")


my_comparisons <- list( c("IDCP", "PCA"))

library(ggpubr)

ggpaired(df_plot_dat, x = "variable", y = "value",shape = "variable",color = "variable",id="hxid",palette = "npg",label="hxid",
         label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Gain Gene CNV numbers")

ggsave(file.path(outdir,paste("sequenza","gain","gene_level_cnv.pdf",sep=".")),width = 6,height = 6)




#### loss

df_plot=df[,c(4,8,22)]
colnames(df_plot)
colnames(df_plot)[1:2]=c("IDCP","PCA")


library(reshape)
df_plot_dat=melt(df_plot,id.vars = "hxid")


my_comparisons <- list( c("IDCP", "PCA"))

library(ggpubr)

ggpaired(df_plot_dat, x = "variable", y = "value",shape = "variable",color = "variable",id="hxid",palette = "npg",label="hxid",
         label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Loss Gene CNV numbers")

ggsave(file.path(outdir,paste("sequenza","loss","gene_level_cnv.pdf",sep=".")),width = 6,height = 6)



ggscatter(df,x="GAIN_CNV_percentage",y="LOSS_CNV_percentage",label = "hxid")+xlim(0,100)+ylim(0,100)
ggsave(file.path(outdir,paste("sequenza","loss_vs_loss","gene_level_cnv.shared.pdf",sep=".")),width = 6,height = 6)



snv=read.csv("02_data_profile/snv.cnv.shared.csv",stringsAsFactors = F,header = T)
snv=snv[,1:7]

## remove hxid 12
snv=snv[snv$hxid!="HXIDCP-12",]


pos=match(snv$hxid,df$hxid)
snv=cbind(snv,df[pos,])

identical(snv[,1],snv[,29])
identical(snv[,2],snv[,26])
identical(snv[,3],snv[,27])

snv[,1:3]=NULL

ggscatter(snv,x="share_snv_percentage",y="TOtal_CNV_percentage",label = "hxid")+xlim(0,100)+ylim(0,100)+xlab("Shared SNV percentage(%)")+ylab("Shared CNV percentage(%)")
ggsave(file.path(outdir,paste("SNV",".sequencza.gene_level_cnv.shared.pdf",sep=".")),width = 6,height = 6)


write.csv(snv[,c("hxid","share_snv_percentage","TOtal_CNV_percentage")],file.path(outdir,"snv_cnv_shared.percentage.csv"),quote = F,row.names = F)


