rm(list=ls())
options(stringsAsFactors = F)

##  setwd("../")

files=list.files(path="10_72425_pyclone",pattern = "loci.tsv",all.files = T,full.names = T,recursive = T)

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
    
    sub_pos=sapply(names(mean_ccf),function(x)which(sub$cluster_id==x))
    
    sub=sub[unlist(sub_pos),]
    
    sub$clonalty="Subclone"
    clone_number=names(mean_ccf)[which.max(mean_ccf)]
    sub$clonalty[sub$cluster_id==clone_number]="Clone"
    
    result=rbind(result,sub)
    
  }
}

result=as.data.frame(result)

result$gene=sapply(result$mutation_id,function(x)unlist(strsplit(x,":"))[5])


samples=read.csv("data/idcp.samples.with.HXid.72425.csv",stringsAsFactors = F)
## 

pos=sapply(samples$Tumor,function(x)which(result$sample_id==x))

result=result[unlist(pos),]

pos=match(result$sample_id,samples$Tumor)
result=cbind(result,samples[pos,])

identical(result$sample_id,result$Tumor)
result$group[result$group=="adenocarcinoma"]="PCA"


result$hxid_group=paste(result$hxid,result$group,sep="_")


write.csv(result,file.path("10_72425_pyclone/","72425.clone_subclone_classification.csv"),quote = F,row.names = F)

