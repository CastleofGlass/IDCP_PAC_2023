rm(list=ls())

#  setwd("../")

dat=read.csv("07_sequenza/sequenza.solutions.munual.csv",stringsAsFactors = F)
dat=dat[dat$rank==1,]

samples=read.csv("data/idcp.samples.with.HXid.update.23samples.csv",stringsAsFactors = F)

pos=match(dat$prefix,samples$Tumor)

dat=cbind(dat,samples[pos,])
dat$group[dat$group=="adenocarcinoma"]="PCA"



length(unique(dat$Tumor))
dat=dat[dat$hxid!="HXIDCP-12",]  # remove hxid-12

length(unique(dat$Tumor))



result_dir="07_sequenza"


library(ggpubr)



my_comparisons <- list( c("IDCP", "PCA"))


  ggpaired(dat, x = "group", y = "cellularity",shape = "group",color = "group",id="hxid",palette = "npg",label="hxid",
           label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
    stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("cellularity")
  ggsave(file.path(result_dir,paste("sequenza","cellularity","pdf",sep=".")),width = 6,height = 6)



write.csv(dat,file.path(result_dir,"sequenza.solutions.update.csv"),quote = F,row.names = F)


ggpaired(dat, x = "group", y = "ploidy",shape = "group",color = "group",id="hxid",palette = "npg",label="hxid",
         label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("ploidy")
ggsave(file.path(result_dir,paste("sequenza","ploidy","pdf",sep=".")),width = 6,height = 6)



library(openxlsx)

origin=openxlsx::read.xlsx("data/snv_cnv_share.xlsx",sheet=1)
origin_b=origin$hxid[origin$Group=="B"]

dat$origin="A_C"

pos=sapply(origin_b,function(x)which(dat$hxid==x))

dat$origin[as.vector(pos)]="B"



dat_shared=dat[dat$origin=="A_C",]
my_comparisons <- list( c("IDCP", "PCA"))

ggpaired(dat_shared, x = "group", y = "ploidy",shape = "group",color = "group",id="hxid",palette = "npg",label="hxid",
         label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("ploidy")
ggsave(file.path(result_dir,paste("sequenza","ploidy.within.A_C","pdf",sep=".")),width = 6,height = 6)



dat_diff=dat[dat$origin=="B",]
my_comparisons <- list( c("IDCP", "PCA"))

ggpaired(dat_diff, x = "group", y = "ploidy",shape = "group",color = "group",id="hxid",palette = "npg",label="hxid",
         label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("ploidy")
ggsave(file.path(result_dir,paste("sequenza","ploidy.within.B","pdf",sep=".")),width = 6,height = 6)




dat_idcp=dat[dat$group=="IDCP",]
my_comparisons <- list( c("A_C", "B"))

ggboxplot(dat_idcp, x = "origin", y = "ploidy",shape = "origin",color = "origin",id="hxid",palette = "npg",label="hxid",
         label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("ploidy")+geom_jitter(shape=16, position=position_jitter(0.2))
ggsave(file.path(result_dir,paste("sequenza","ploidy.within.IDCP","pdf",sep=".")),width = 6,height = 6)




dat_pca=dat[dat$group=="PCA",]
my_comparisons <- list( c("A_C", "B"))

ggboxplot(dat_pca, x = "origin", y = "ploidy",shape = "origin",color = "origin",id="hxid",palette = "npg",label="hxid",
          label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("ploidy")+geom_jitter(shape=16, position=position_jitter(0.2))
ggsave(file.path(result_dir,paste("sequenza","ploidy.within.PCA","pdf",sep=".")),width = 6,height = 6)









