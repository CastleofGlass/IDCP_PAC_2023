rm(list=ls())
options(stringsAsFactors = F)

##  setwd("../")

files=list.files(path="05_run_pyclone",pattern = "loci.tsv",all.files = T,full.names = T,recursive = T)

result=c()


for(file in files){
  dat=read.table(file,stringsAsFactors = F,header = T,sep="\t")
  sample_ids=unique(dat$sample_id)

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


library(ggpubr)

my_comparisons <- list( c("IDCP", "PCA"))
ggboxplot(result, x = "group", y = "Cluster_Num",shape = "group",color = "group",notch = FALSE,palette = "jco",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+ylab("Cluster Number")

ggsave(file.path("05_run_pyclone","05.pyclone.clone.number.pdf"),width = 6,height = 6)




