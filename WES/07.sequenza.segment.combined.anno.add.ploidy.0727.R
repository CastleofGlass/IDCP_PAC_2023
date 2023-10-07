rm(list=ls())

# setwd("../")

outdir="07_sequenza"


cnv=read.csv("07_sequenza/sequenza.segment.combined.anno.csv",stringsAsFactors = F,check.names = F)

solutions=read.csv("07_sequenza/sequenza.solutions.munual.csv",stringsAsFactors = F,check.names = F)
solutions=solutions[solutions$rank==1,]

pos=match(cnv$tumor,solutions$prefix)

cnv=cbind(cnv,solutions[pos,c(1,2,4)])



filter_files=list.files(path="07_sequenza/cnv_anno/",pattern = ".sequenza.result.anno.bed.filtered.tsv",all.files = T,full.names = T)

filter_df=c()

for(file in filter_files){
  cnv_prefix=gsub(".sequenza.result.anno.bed.filtered.tsv","",basename(file))
  dat=read.table(file,stringsAsFactors = F,check.names = F,sep="\t",header = T)
  dat$cnv_prefix=cnv_prefix
  
  filter_df=rbind(filter_df,dat)
  
}


filter_df$newid=paste(filter_df$chromosome,filter_df$start.pos,filter_df$end.pos,filter_df$CNt,filter_df$cnv_prefix,sep="_")


cnv$match=paste(cnv$chromosome,cnv$start.pos,cnv$end.pos,cnv$CNt,cnv$prefix,sep="_")

newpos=match(cnv$match,filter_df$newid)


cnv=cbind(cnv,filter_df[newpos,c("CNV_classification","cnv_prefix","newid")])


cnv$cnv_prefix=NULL
cnv$newid=NULL
cnv$match=NULL


write.csv(cnv,file.path(outdir,"sequenza.segment.combined.anno.add.ploidy.20210727.csv"),quote = F,row.names = F)







