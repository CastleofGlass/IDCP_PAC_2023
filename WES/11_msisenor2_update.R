#  setwd("../data/")

dat=read.csv("msisensor.summary.csv",stringsAsFactors = F)

samples=read.csv("idcp.samples.with.HXid.update.23samples.csv",stringsAsFactors = F)


samples$group[samples$group=="adenocarcinoma"]="PCA"


pos=match(dat$Sample,samples$Tumor)

res=cbind(dat,samples[pos,])

write.csv(res,"msisensor.result.csv",quote = F,row.names = F)
