##setwd("../")

rm(list=ls())

library(magrittr)
library(openxlsx)

dat=openxlsx::read.xlsx("data/20210413_checked_21h.xlsx", sheet=1)
colnames(dat)[c(1,12,13)]=c("CP","transcript","variant_class")

dim(dat)  #6713
table(dat$variant_class)
# CNV 3131
# SNV  3582

result_dir="01_data_preprocess"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}




### for CNV:
cnv=dat[dat$variant_class=="CNV",]
cnv_used=cnv[,c("CP","sample_id","gene","chr","chr_start","chr_end","freq","p_dot")]
colnames(cnv_used)[7:8]=c("cn","cnv_type")

write.csv(cnv_used,file.path("data","0413.cnv.csv"),quote = F,row.names = F)



samples=unique(cnv_used$sample_id)  #67 samples

cnv_number=c()

for(sample in samples){
  
  outfile=file.path(result_dir,paste(sample,".cnv.tsv",sep=""))
  sub=cnv_used[cnv_used$sample_id==sample,]
  
  write.table(sub,outfile,quote = F,sep="\t",row.names = F)
  
  cnv_number=rbind(cnv_number,c(sample,nrow(sub),nrow(sub[sub$cnv_type=="loss",]), nrow(sub[sub$cnv_type=="gain",])))
}

cnv_number=as.data.frame(cnv_number)
colnames(cnv_number)=c("Sample_id","Number_of_total_CNV","Number_of_CNV_loss","Number_of_CNV_gain")


write.csv(cnv_number,file.path("data","0413.cnv.number.csv"),quote = F,row.names = F)




###########for SNV

snvs=dat[dat$variant_class=="SNV",]  ##3582 somatics SNVs
write.csv(snvs,file.path("data","0413.snv.csv"),quote = F,row.names = F)

library(ggpubr)
p=ggdensity(snvs,x="freq",rug = TRUE)
ggsave(file.path(result_dir,"01_all_80_snv_vaf_distribution.pdf"))

samples=unique(snvs$sample_id)  #80 samples

snv_number=c()

for(sample in samples){
  
  outfile=file.path(result_dir,paste(sample,".snv.tsv",sep=""))
  sub=snvs[snvs$sample_id==sample,]
  
  write.table(sub,outfile,quote = F,sep="\t",row.names = F)
  
  snv_number=rbind(snv_number,c(sample,nrow(sub)))
}

snv_number=as.data.frame(snv_number)
colnames(snv_number)=c("Sample_id","Number_of_total_SNV")


write.csv(snv_number,file.path("data","0413.snv.number.csv"),quote = F,row.names = F)


