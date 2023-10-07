rm(list=ls())
options(stringsAsFactors = F)

##  setwd("../")

files=list.files(path="09_pyclone",pattern = "loci.tsv",all.files = T,full.names = T,recursive = T)

result=c()


for(file in files){
  dat=read.table(file,stringsAsFactors = F,header = T,sep="\t")
  sample_ids=unique(dat$sample_id)
  
  cluster_vars=table(dat$cluster_id)
  clusters=cluster_vars[cluster_vars>=3*length(sample_ids)]  ## at least 3 variants in the clusters
  
  pos=sapply(names(clusters),function(x)which(dat$cluster_id==x))
  
  print(nrow(dat))
  dat=dat[unlist(pos),]
  print(nrow(dat))
  
  
  for(id in sample_ids){
    
    sub=dat[dat$sample_id==id,]
    
    mean_ccf=tapply(sub$cellular_prevalence,sub$cluster_id,mean)
    mean_ccf=mean_ccf[mean_ccf>0.05]
    result=rbind(result,c(id,length(mean_ccf)))
  }
}

result=as.data.frame(result)
colnames(result)=c("ID","Cluster_Num")

anno=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F)

pos=match(result$ID,anno$Tumor)
result=cbind(result,anno[pos,])

result$Cluster_Num=as.numeric(result$Cluster_Num)

result$group[result$group=="adenocarcinoma"]="PCA"

## remove hxid-12
result=result[!is.na(result$hxid),]

dim(result)

write.csv(result,file.path("09_pyclone","09.pyclone.clone.number.csv"),quote = F,row.names = F)

library(ggpubr)

my_comparisons <- list( c("IDCP", "PCA"))

result$xj <- jitter(result$Cluster_Num)

ggpaired(result, x = "group", y = "xj",shape = "group",color = "group",id="hxid",palette = "npg",label="hxid",
         label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash",ylab="Cluster Number")+
  #geom_jitter(width = 0.1,height = 0.1)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)

ggsave(file.path("09_pyclone","09.pyclone.clone.number.pdf"),width = 6,height = 6)




