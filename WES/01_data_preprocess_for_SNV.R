## setwd("../")

rm(list=ls())

library(magrittr)
library(openxlsx)
library(dplyr)

dat=openxlsx::read.xlsx("data/20210413_checked_21h.xlsx", sheet=1)
colnames(dat)[c(1,12,13)]=c("CP","transcript","variant_class")

dim(dat)  #6713
table(dat$variant_class)

result_dir="01_data_preprocess"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}


snvs=dat[dat$variant_class=="SNV",]  ##3582 somatics SNVs

samples=read.csv("data/idcp.samples.with.HXid.update.23samples.csv",stringsAsFactors = F)

samples_used=intersect(samples$Tumor,snvs$sample_id)


snvs_res <- snvs %>%
  dplyr::filter(sample_id %in% samples_used)


pos=match(snvs_res$sample_id,samples$Tumor)
snvs_result=cbind(snvs_res,samples[pos,])

write.csv(snvs_result,file.path(result_dir,"20210413_snvs_with_hxid.csv"),quote = F,row.names = F)

