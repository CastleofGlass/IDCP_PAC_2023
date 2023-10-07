rm(list=ls())

# setwd("../")
library(SomaticSignatures)
library(SomaticCancerAlterations)
#library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)


### read snvs
snvs=read.csv("data/0413.snv.csv",stringsAsFactors = F)

groups=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F)
groups$group[groups$group=="adenocarcinoma"]="PCA"

groups$hxid_group=paste(groups$hxid,groups$group,sep="_")

pos=sapply(groups$Tumor,function(x)which(snvs$sample_id==x))

dim(snvs)
snvs=snvs[unlist(pos),]
dim(snvs)

snvs=snvs[nchar(snvs$reference_allele)==1 & nchar(snvs$variant_allele)==1,]
snvs=snvs[snvs$reference_allele!="-",]
snvs=snvs[snvs$variant_allele!="-",]


dim(snvs)


pos=match(snvs$sample_id,groups$Tumor)
snvs=cbind(snvs,groups[pos,])

identical(snvs$sample_id,snvs$Tumor)


result_dir="08_mutation_signatures"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}


sca_vr = VRanges(
  seqnames =  snvs$chr ,
  ranges = IRanges(start = snvs$chr_start,end = snvs$chr_end),
  ref = snvs$reference_allele,
  alt = snvs$variant_allele,
  sampleNames = snvs$hxid_group,
  study=snvs$group
  )

sca_vr


sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19)


##
sca_mm=motifMatrix(sca_motifs, group = "sampleNames", normalize = FALSE)


rename_fun=function(x){
 ref=substr(x,1,1)
 alt=substr(x,2,2)
 up=substr(x,4,4)
 down=substr(x,6,6)
 
 res=paste(up,"[",ref,">",alt,"]",down,sep="")
 return(res)
}


sca_mm_rowname=sapply(rownames(sca_mm),rename_fun)

df=cbind(sca_mm_rowname,sca_mm)
colnames(df)[1]="Mutation Types"


write.table(df,file.path("08_mutation_signatures","Matrix.txt"),sep="\t",quote = F,row.names = F)



