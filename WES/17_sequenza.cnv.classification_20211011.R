rm(list=ls())

#  setwd("../")

dat=read.csv("07_sequenza/sequenza.segment.combined.anno.add.ploidy.20210727.csv",stringsAsFactors = F)
dat$group[dat$group=="adenocarcinoma"]="PCA"

length(unique(dat$Tumor))
dat=dat[dat$hxid!="HXIDCP-12",]  # remove hxid-12

length(unique(dat$Tumor))



result_dir="17_sequenza_logR_0.35"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}



dat$logR=log2(dat$depth.ratio)

dat$logR_CNV_classification="."
dat$logR_CNV_classification[dat$logR>0.35]="gain"
dat$logR_CNV_classification[dat$logR< (-0.35)]="loss"


table(dat$CNV_classification,dat$logR_CNV_classification)

write.csv(dat,file.path(result_dir,"sequenza.segment.combined.anno.add.ploidy.20211011.csv"),quote = F,row.names = F)





############# AR PTEN check:
sequenza=dat

contain_AR=sapply(sequenza$Genes,function(x){y=unlist(strsplit(x,";"));z=length(intersect("AR",y))==1;return(z)})

sequenza_AR=sequenza[contain_AR,]
sequenza_AR$NOTE=paste(sequenza_AR$logR_CNV_classification," (",sequenza_AR$CNt,",", round(sequenza_AR$logR,3),") ",sep="")


contain_pten=sapply(sequenza$Genes,function(x){y=unlist(strsplit(x,";"));z=length(intersect("PTEN",y))==1;return(z)})

sequenza_PTEN=sequenza[contain_pten,]
sequenza_PTEN$NOTE=paste(sequenza_PTEN$logR_CNV_classification," (",sequenza_PTEN$CNt,",", round(sequenza_PTEN$logR,3),") ",sep="")



report=read.csv("data/0413.cnv.csv",stringsAsFactors = F,header = T)

report$NOTE=paste(report$cnv_type," (",report$cn,") ",sep="")

##### filter genes

AR=report[report$gene=="AR",]
PTEN=report[report$gene=="PTEN",]



samples=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F)
samples$group[samples$group=="adenocarcinoma"]="PCA"

AR_pos=match(samples$Tumor,AR$sample_id)
samples$AR_report=AR$NOTE[AR_pos]

PTEN_pos=match(samples$Tumor,PTEN$sample_id)
samples$PTEN_report=PTEN$NOTE[PTEN_pos]



#### 
sequenza_AR_pos=match(samples$Tumor,sequenza_AR$prefix)
samples$AR_sequenza=sequenza_AR$NOTE[sequenza_AR_pos]


sequenza_PTEN_pos=match(samples$Tumor,sequenza_PTEN$prefix)
samples$PTEN_sequenza=sequenza_PTEN$NOTE[sequenza_PTEN_pos]


pos=match(samples$Tumor,sequenza$prefix)
samples$ploidy_sequenza=sequenza$ploidy[pos]


library("openxlsx")

openxlsx::write.xlsx(samples,file.path(result_dir,"AR.PTEN.check.result.20211011.xlsx"),overwrite = TRUE)




rm(list=ls())
options(stringsAsFactors = F)

#setwd("../")

result_dir="17_sequenza_logR_0.35"

samples=read.csv("data/idcp.pairs.with.HXid.csv",stringsAsFactors = F)

cnv=read.csv(file.path(result_dir,"sequenza.segment.combined.anno.add.ploidy.20211011.csv"),stringsAsFactors = F)



result=c()


myfun=function(df){
  dat=c()
  for(n in 1:nrow(df)){
    genes=df$Genes[n]
    genes=unlist(strsplit(genes,";"))
    cnv_type=df$logR_CNV_classification[n]
    
    tmp=paste(genes,cnv_type,sep="@")
    dat=c(dat,tmp)
  }
  return(dat)
}


for(i in 1:nrow(samples)){
  print(i)
  idcp_id=samples$idcp[i]
  ade_id=samples$adenocarcinoma[i]
  
  idcp_cnv=cnv[cnv$prefix==idcp_id,]
  idcp_cnv=idcp_cnv[idcp_cnv$logR_CNV_classification!=".",]
  
  ade_cnv=cnv[cnv$prefix==ade_id,]
  ade_cnv=ade_cnv[ade_cnv$logR_CNV_classification!=".",]
  
  
  idcp_segment_gain_number=nrow(idcp_cnv[idcp_cnv$logR_CNV_classification=="gain",])
  idcp_segment_loss_number=nrow(idcp_cnv[idcp_cnv$logR_CNV_classification=="loss",])
  
  
  ade_segment_gain_number=nrow(ade_cnv[ade_cnv$logR_CNV_classification=="gain",])
  ade_segment_loss_number=nrow(ade_cnv[ade_cnv$logR_CNV_classification=="loss",])
  
  
  idcp_gene_classification=myfun(idcp_cnv)
  pca_gene_classification=myfun(ade_cnv)
  
  
  total_cnv=unique(c(idcp_gene_classification,pca_gene_classification))
  
  share_cnv=intersect(idcp_gene_classification,pca_gene_classification)
  
  share_cnv_gain=share_cnv[grep("@gain",share_cnv)]
  share_cnv_loss=share_cnv[grep("@loss",share_cnv)]
  
  res=c(idcp_id,nrow(idcp_cnv),idcp_segment_gain_number,idcp_segment_loss_number,ade_id,nrow(ade_cnv),ade_segment_gain_number,ade_segment_loss_number,
        length(total_cnv),length(share_cnv),length(share_cnv_gain),length(share_cnv_loss))
  
  result=rbind(result,res)
}





result=as.data.frame(result)
colnames(result)=c("IDCP_ID","IDCP_Number_of_CNV_Segment","IDCP_Number_of_GAIN_Segment","IDCP_Number_of_LOSS_Segment",
                   "PCA_ID","PCA_Number_of_CNV_Segment","PCA_Number_of_GAIN_Segment","PCA_Number_of_LOSS_Segment",
                   "Total_CNV_gene_level","Total_CNV_Overlapped_gene_level",
                   "Gain_CNV_Overlapped_gene_level","Loss_CNV_overlapped_level")



outdir=result_dir


write.csv(result,file.path(outdir,"sequenza_gene_level_cnv_overlapped.csv"),quote = F,row.names = F)



df=read.csv(file.path(outdir,"sequenza_gene_level_cnv_overlapped.csv"),stringsAsFactors = F,header = T)

df$TOtal_CNV_percentage=100*df$Total_CNV_Overlapped_gene_level/df$Total_CNV_gene_level


pos=match(df$IDCP_ID,samples$idcp)
df=cbind(df,samples[pos,])


write.csv(df,file.path(outdir,"sequenza_gene_level_cnv_overlapped_update.csv"),quote = F,row.names = F)



##################SNV shared
snv=read.csv("02_data_profile/snv.cnv.shared.csv",stringsAsFactors = F,header = T)
snv=snv[,1:7]

## remove hxid 12
snv=snv[snv$hxid!="HXIDCP-12",]


pos=match(snv$hxid,df$hxid)
snv=cbind(snv,df[pos,])

identical(snv[,1],snv[,25])
identical(snv[,2],snv[,22])
identical(snv[,3],snv[,23])



df=snv[,c("share_snv_percentage","TOtal_CNV_percentage")]
rownames(df)=snv$hxid

pdf(file.path(outdir,"snv.cnv.shared.hca.20211011.pdf",sep=""),width = 6,height = 6)
plot(hclust(dist(df,method = "euclidean"),method = "ward.D"))
dev.off()


library(ggpubr)

snv_df=df
snv_df$ID=rownames(snv_df)
snv_df$newid=gsub("HXIDCP-","",snv_df$ID)

ggscatter(snv_df,x="share_snv_percentage",y="TOtal_CNV_percentage",label = "newid",xlab="Shared SNV percentage (%)",
          ylab="Shared CNV percentage (%)")

ggsave(file.path(outdir,"snv.cnv.shared.scatter.plot.20211011.pdf"),width = 7,height = 7)






