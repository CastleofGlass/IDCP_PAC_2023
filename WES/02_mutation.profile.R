rm(list=ls())

setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ/")

samples=read.table("data/80.samples.list",stringsAsFactors = F,header = T)

hxid=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F,header = T)  ##46 samples
hxid$group[hxid$group=="adenocarcinoma"]="PCA"

###  variants number check:

snv_num=read.csv("data/0413.snv.number.csv",stringsAsFactors = F)
cnv_num=read.csv("data/0413.cnv.number.csv",stringsAsFactors = F)

pos_snv=match(hxid$Tumor,snv_num$Sample_id)
pos_cnv=match(hxid$Tumor,cnv_num$Sample_id)

dat=cbind(hxid,snv_num[pos_snv,],cnv_num[pos_cnv,])

identical(dat[,2],dat[,6])
table(dat[,2]==dat[,8])


colnames(dat)[8]="sample_id_2"

dat$group[dat$group=="adenocarcinoma"]="PCA"



result_dir="02_data_profile"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}






library(ggpubr)

### FOR SNV
p=ggbarplot(dat,x="hxid","Number_of_total_SNV",fill="group",color = "group",palette = "Paired",
            position = position_dodge(0.9),x.text.angle=90,rotate=TRUE,ggtheme = theme_minimal(),label=FALSE
)

p+ylab("Number of Mutations")+xlab("")
ggsave(file.path(result_dir,"02.snv.number.pdf"),width = 8,height = 8)



### for CNV

p=ggbarplot(dat,x="hxid","Number_of_total_CNV",fill="group",color = "group",palette = "Paired",
            position = position_dodge(0.9),x.text.angle=90,rotate=TRUE,ggtheme = theme_minimal(),label=FALSE
)

p+ylab("Number of Total CNV")+xlab("")
ggsave(file.path(result_dir,"02.CNV.number.pdf"),width = 8,height = 8)



### for CNV GAIN

p=ggbarplot(dat,x="hxid","Number_of_CNV_gain",fill="group",color = "group",palette = "Paired",
            position = position_dodge(0.9),x.text.angle=90,rotate=TRUE,ggtheme = theme_minimal(),label=FALSE
)

p+ylab("Number of CNV Gain")+xlab("")
ggsave(file.path(result_dir,"02.CNV.number.gain.pdf"),width = 8,height = 8)


### for CNV loss

p=ggbarplot(dat,x="hxid","Number_of_CNV_loss",fill="group",color = "group",palette = "Paired",
            position = position_dodge(0.9),x.text.angle=90,rotate=TRUE,ggtheme = theme_minimal(),label=FALSE
)

p+ylab("Number of CNV Loss")+xlab("")
ggsave(file.path(result_dir,"02.CNV.number.loss.pdf"),width = 8,height = 8)





##### idcp vs pca
my_comparisons <- list( c("IDCP", "PCA"))
ggboxplot(dat, x = "group", y = "Number_of_total_SNV",
          add = "jitter", shape = "group",color = "group",notch = FALSE)+
          stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Number of Total SNV")

ggsave(file.path(result_dir,"02.snv.idcp_vs_pca.pdf"),width = 6,height = 6)


ggboxplot(dat, x = "group", y = "Number_of_total_CNV",
          add = "jitter", shape = "group",color = "group",notch = FALSE)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Number of Total CNV")

ggsave(file.path(result_dir,"02.CNV.idcp_vs_pca.pdf"),width = 6,height = 6)



####
snvs=read.csv("data/0413.snv.csv",stringsAsFactors = F,header = T)
snvs$variant_id=paste(snvs$chr,snvs$chr_start,snvs$chr_end,snvs$reference_allele,snvs$variant_allele,sep="_")
head(snvs)

cnvs=read.csv("data/0413.cnv.csv",stringsAsFactors = F,header = T)
cnvs$variant_id=paste(cnvs$gene,cnvs$cnv_type,sep="_")
head(cnvs)

ids=unique(hxid$hxid)
ids

head(hxid)

variant_share=c()

for(id in ids){
  idcp_id=hxid$Tumor[hxid$hxid==id&hxid$group=="IDCP"]
  pca_id=hxid$Tumor[hxid$hxid==id&hxid$group=="PCA"]
 
  idcp_snv=snvs[snvs$sample_id==idcp_id,]
  pca_snv=snvs[snvs$sample_id==pca_id,]
  
  idcp_cnv=cnvs[cnvs$sample_id==idcp_id,]
  pca_cnv=cnvs[cnvs$sample_id==pca_id,]
  
  res=c(id,idcp_id,pca_id,nrow(idcp_snv),nrow(pca_snv),length(intersect(idcp_snv$variant_id,pca_snv$variant_id)),100*length(intersect(idcp_snv$variant_id,pca_snv$variant_id))/length(union(idcp_snv$variant_id,pca_snv$variant_id)),
                          nrow(idcp_cnv),nrow(pca_cnv),length(intersect(idcp_cnv$variant_id,pca_cnv$variant_id)),100*length(intersect(idcp_cnv$variant_id,pca_cnv$variant_id))/length(union(idcp_cnv$variant_id,pca_cnv$variant_id)))
  variant_share=rbind(variant_share,res)
}

variant_share=as.data.frame(variant_share)
colnames(variant_share)=c("hxid","idcp","pca","IDCP_snv_number","pca_snv_number","share_snv_number","share_snv_percentage","IDCP_cnv_number","pca_cnv_number","share_cnv_number","share_cnv_percentage")
head(variant_share)

variant_share[,4:ncol(variant_share)]=apply(variant_share[,4:ncol(variant_share)],2,as.numeric)
head(variant_share)

write.csv(variant_share,file.path("02_data_profile","snv.cnv.shared.csv"),quote = F,row.names = F)


ggscatter(variant_share,x="share_snv_percentage",y="share_cnv_percentage",xlab="SNV shared Percentage",ylab="CNV shared percentage",
          label="hxid")

ggsave(file.path(result_dir,"02.snv_vs_cnv.shared.percentage.pdf"),width = 9,height = 9)







