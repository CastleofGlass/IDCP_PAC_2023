#  setwd("../data")

rm(list=ls())

dat=read.table("scarHRD.summary.tsv",stringsAsFactors = F,sep="\t",header = T,check.names = F)

groups=read.csv("idcp.samples.with.HXid.update.23samples.csv",stringsAsFactors = F)
groups$group[groups$group=="adenocarcinoma"]="PCA"


pos=match(dat$tumor,groups$Tumor)


dat=cbind(dat,groups[pos,])

if(identical(dat$tumor,dat$Tumor)){
  print("duplicated tumor colname")
  dat$Tumor=NULL
}

write.csv(dat,"scarHRD.summary.csv",quote = F,row.names = F)






samples=read.csv("idcp.samples.with.HXid.update.23samples.csv",stringsAsFactors = F)


samples$group[samples$group=="adenocarcinoma"]="PCA"


pos=match(dat$Sample,samples$Tumor)

res=cbind(dat,samples[pos,])

write.csv(res,"msisensor.result.csv",quote = F,row.names = F)
