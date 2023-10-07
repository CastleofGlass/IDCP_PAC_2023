#  setwd("../")


rm(list=ls())


library(openxlsx)


newgroup=openxlsx::read.xlsx("data/IDCP.groups.by.sun.xlsx",sheet=1)
newgroup=newgroup[,c(1,2,5)]
colnames(newgroup)=c("HXID","ADT","Mutation_group")

samples=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F)


pos=match(samples$hxid,newgroup$HXID)

group=cbind(samples,newgroup[pos,])


outdir="12_gistic"

if(!dir.exists(outdir)){
  dir.create(outdir)
}

write.csv(group,file.path(outdir,"12.gistic.sample.group.csv"),quote = F,row.names = F)




