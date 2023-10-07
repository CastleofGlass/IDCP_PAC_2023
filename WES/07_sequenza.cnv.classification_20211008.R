rm(list=ls())

#  setwd("../")

dat=read.csv("07_sequenza/sequenza.segment.combined.anno.add.ploidy.20210727.csv",stringsAsFactors = F)
dat$group[dat$group=="adenocarcinoma"]="PCA"

length(unique(dat$Tumor))
dat=dat[dat$hxid!="HXIDCP-12",]  # remove hxid-12

length(unique(dat$Tumor))



result_dir="07_sequenza"

dat$logR=log2(dat$depth.ratio)

dat$logR_CNV_classification="."
dat$logR_CNV_classification[dat$logR>log2(1.5)]="gain"
dat$logR_CNV_classification[dat$logR< (-log2(1.5))]="loss"


table(dat$CNV_classification,dat$logR_CNV_classification)

write.csv(dat,"07_sequenza/sequenza.segment.combined.anno.add.ploidy.20211009.csv",quote = F,row.names = F)





############# AR PTEN check:
sequenza=dat

contain_AR=sapply(sequenza$Genes,function(x){y=unlist(strsplit(x,";"));z=length(intersect("AR",y))==1;return(z)})

sequenza_AR=sequenza[contain_AR,]
sequenza_AR$NOTE=paste(sequenza_AR$logR_CNV_classification," (",sequenza_AR$CNt,",", round(sequenza_AR$logR,3),") ",sep="")


contain_pten=sapply(sequenza$Genes,function(x){y=unlist(strsplit(x,";"));z=length(intersect("PTEN",y))==1;return(z)})

sequenza_PTEN=sequenza[contain_pten,]
sequenza_PTEN$NOTE=paste(sequenza_PTEN$logR_CNV_classification," (",sequenza_PTEN$CNt,",", round(sequenza_PTEN$logR,3),") ",sep="")



report=read.csv("data/0413.cnv.csv",stringsAsFactors = F,header = T)

report$NOTE=paste(report$cnv_type," (",report$cn,") ",sep="")

##### filter genes

AR=report[report$gene=="AR",]
PTEN=report[report$gene=="PTEN",]



samples=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F)
samples$group[samples$group=="adenocarcinoma"]="PCA"

AR_pos=match(samples$Tumor,AR$sample_id)
samples$AR_report=AR$NOTE[AR_pos]

PTEN_pos=match(samples$Tumor,PTEN$sample_id)
samples$PTEN_report=PTEN$NOTE[PTEN_pos]



#### 
sequenza_AR_pos=match(samples$Tumor,sequenza_AR$prefix)
samples$AR_sequenza=sequenza_AR$NOTE[sequenza_AR_pos]


sequenza_PTEN_pos=match(samples$Tumor,sequenza_PTEN$prefix)
samples$PTEN_sequenza=sequenza_PTEN$NOTE[sequenza_PTEN_pos]


pos=match(samples$Tumor,sequenza$prefix)
samples$ploidy_sequenza=sequenza$ploidy[pos]


library("openxlsx")

openxlsx::write.xlsx(samples,"AR.PTEN.check.result.20211009.xlsx",overwrite = TRUE)




