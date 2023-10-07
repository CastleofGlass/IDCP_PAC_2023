rm(list=ls())

#  setwd("../")

dat=read.csv("07_sequenza/sequenza.segment.combined.csv",stringsAsFactors = F)
dat$group[dat$group=="adenocarcinoma"]="PCA"

length(unique(dat$Tumor))
dat=dat[dat$hxid!="HXIDCP-12",]  # remove hxid-12

length(unique(dat$Tumor))



result_dir="07_sequenza"


library(ggpubr)

ggscatter(dat,x="N.BAF","N.ratio",color ="tumor" )

ggsave(file.path(result_dir,"07.N.BAF_vs_N.ratio.pdf"),width = 8,height = 10)



ggscatter(dat,x="N.BAF","Bf",color ="tumor" )

ggsave(file.path(result_dir,"07.N.BAF_vs_Bf.pdf"),width = 8,height = 10)

ggscatter(dat,x="N.ratio","depth.ratio",color ="tumor" )

ggsave(file.path(result_dir,"07.N.ratio_vs_depth.ratio.pdf"),width = 8,height = 10)



ggdensity(dat,x="N.BAF",rug = TRUE)+scale_x_log10()
ggsave(file.path(result_dir,"07.N.BAF.density.pdf"),width = 7,height = 7)


ggdensity(dat,x="N.ratio",rug = TRUE)+scale_x_log10()
ggsave(file.path(result_dir,"07.N.ratio.density.pdf"),width = 7,height = 7)

res=dat[dat$N.BAF>=9&dat$N.ratio>=30,]

res$HX_group=paste(res$hxid,res$group,sep="_")

res$cnv_type="Normal"
res$cnv_type[res$CNt>2]="Gain"
res$cnv_type[res$CNt<2]="Loss"


ids=unique(res$HX_group)

ids_df=c()

for(id in ids){
  sub=res[res$HX_group==id,]
  sub=sub[sub$cnv_type!="Normal",]
  
  all=sum(sub$end.pos-sub$start.pos)/1000000
  
  gain=sub[sub$cnv_type=="Gain",]
  
  gain_num=sum(gain$end.pos-gain$start.pos)/1000000
  
  loss=sub[sub$cnv_type=="Loss",]
  loss_num=sum(loss$end.pos-loss$start.pos)/1000000
  
  ids_df=rbind(ids_df,c(id,nrow(sub),all,nrow(gain),gain_num,nrow(loss),loss_num))
}


ids_df=as.data.frame(ids_df,stringsAsFactors = F)
colnames(ids_df)=c("HXID_group","CNV_number","CNV_burden","Gain_number","Gain_burden","Loss_number","Loss_burden")

ids_df[2:ncol(ids_df)]=apply(ids_df[2:ncol(ids_df)],2,as.numeric)


ids_df$hxid=sapply(ids_df$HXID_group,function(x)unlist(strsplit(x,"_"))[1])
ids_df$group=sapply(ids_df$HXID_group,function(x)unlist(strsplit(x,"_"))[2])




my_comparisons <- list( c("IDCP", "PCA"))

for(variant in colnames(ids_df)[2:7]){
  ggpaired(ids_df, x = "group", y = variant,shape = "group",color = "group",id="hxid",palette = "npg",label="hxid",
           label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
    stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab(variant)
  ggsave(file.path(result_dir,paste("sequenza",variant,"pdf",sep=".")),width = 6,height = 6)
}


write.csv(ids_df,file.path(result_dir,"sequenza.cnv.summary.csv"),quote = F,row.names = F)

