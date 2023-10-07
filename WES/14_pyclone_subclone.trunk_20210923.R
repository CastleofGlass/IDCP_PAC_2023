rm(list=ls())
options(stringsAsFactors = F)

##  setwd("../")
result_dir="09_pyclone"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}



snvs=read.csv("09_pyclone/clone_subclone_classification.csv",stringsAsFactors = F)

drivers=read.table("mobster_378_driver_genes_from_Martincorena_and_Tarabichi.txt",stringsAsFactors = F)

is.driver=function(x){
  y=intersect(x,drivers$V1)
  
  if(length(y)>0){
    z="TRUE"
  }else{
    z="FALSE"
  }
  return(z)
  
}


snvs$is.driver.gene=sapply(snvs$gene,is.driver)

snvs$Trunk_or_Branch="Branch"


snvs$Trunk_or_Branch[snvs$is.driver.gene=="FALSE"]=""

write.csv(snvs,file.path(result_dir,"clone_subclone_classification.add.driver.20210923.csv"),quote = F,row.names = F)