

# Date: 2021-12-13
# Removed CNVs located on chrX, chrY and chrM
# 宋博原来的代码中的nrow 存在问题，会产生NA，倒是number of CNV增加
# 注意不要使用Tumor ID，使用HXID

rm(list=ls())

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyr)
library(agricolae)
library(multcomp)

setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ")

outdir <- "07_sequenza.CNV.burden_20211214"
if(!dir.exists(outdir)){
  dir.create(outdir)
}

chr_len <- read.table("scripts/hg19.fai.txt")
chr_len
colnames(chr_len) <- c("Chromosome", "length")


# 44+2 淋巴转移
dat=read.csv("07_sequenza/sequenza.segment.combined.anno.add.ploidy.20210727.csv",stringsAsFactors = F)
dat[1:2,]
dim(dat)
## [1] 15239    24

dat$group[dat$group == "adenocarcinoma"] <- "PCA"
dat$full_hxid <- paste(dat$hxid, dat$group, sep = "_")

# 不看chrX, chrY, chrM 
dat <- dat %>%  dplyr::filter(chromosome != "chrM" & chromosome != "chrX" & chromosome != "chrY")
dim(dat)
## [1] 14904    24

table(dat$chromosome)

chromosomes <- unique(dat$chromosome) 
chromosomes
samples=unique(dat$full_hxid)
samples





############################### Part 1 ############################### 
# 分染色体
cnv_burden_chr <- c()
# Total CNV length
for (sample in samples) {
  for (chromosome in chromosomes) {
    sub_chr_gain <- dat[which(dat$chromosome == chromosome & dat$full_hxid == sample & dat$CNV_classification == "gain"), ] 
    sub_chr_gain_length <- sum(sub_chr_gain$end.pos - sub_chr_gain$start.pos + 1)
    sub_chr_gain_prop <- round(sub_chr_gain_length / chr_len[which(chr_len$Chromosome == chromosome),"length"], 5)
    
    sub_chr_loss <- dat[which(dat$chromosome == chromosome & dat$full_hxid == sample & dat$CNV_classification == "loss"), ] 
    sub_chr_loss_length <- sum(sub_chr_loss$end.pos - sub_chr_loss$start.pos + 1)
    sub_chr_loss_prop <- round(sub_chr_loss_length / chr_len[which(chr_len$Chromosome == chromosome),"length"], 5)
    
    sub_chr_CNV_length <- sub_chr_loss_length + sub_chr_gain_length
    sub_chr_CNV_prop <- round(sub_chr_CNV_length / chr_len[which(chr_len$Chromosome == chromosome),"length"],5)
    
    cnv_burden_chr <- rbind(cnv_burden_chr, c(sample, chromosome, sub_chr_gain_length, sub_chr_loss_length, sub_chr_CNV_length, sub_chr_gain_prop, sub_chr_loss_prop, sub_chr_CNV_prop))
  }
}

cnv_burden_chr <- as.data.frame(cnv_burden_chr)
colnames(cnv_burden_chr) <- c("HXID", "Chr", "Gian_length", "Loss_length", "CNV_length", "Gian_length_prop", "Loss_length_prop", "CNV_length_prop")
str(cnv_burden_chr)
dim(cnv_burden_chr)
## [1] 1012    8
cnv_burden_chr$Chr <- as.factor(cnv_burden_chr$Chr) # 方差分析需要是因子
cnv_burden_chr$Gian_length <- as.numeric(cnv_burden_chr$Gian_length)
cnv_burden_chr$Loss_length <- as.numeric(cnv_burden_chr$Loss_length)
cnv_burden_chr$Gian_length_prop <- as.numeric(cnv_burden_chr$Gian_length_prop)
cnv_burden_chr$Loss_length_prop <- as.numeric(cnv_burden_chr$Loss_length_prop)
cnv_burden_chr$CNV_length_prop <- as.numeric(cnv_burden_chr$CNV_length_prop)
cnv_burden_chr[1:10,]


write.csv(cnv_burden_chr, file.path(outdir, "CNV_burden_20211214_HXID.csv"), quote = F, row.names = F)

library(ggpubr)
ggboxplot(cnv_burden_chr, x = "Chr", y = "Loss_length_prop")
summary(cnv_burden_chr)


################
## 1. CNV burden: CNV gain + loss 
{
  # 使用spread函数将gd1_long长数据转换为宽数据gd1_wide
  cnv_burden_chr_CNV_prop_broad <- spread(cnv_burden_chr[, c("HXID", "Chr", "CNV_length_prop")], Chr, CNV_length_prop)
  head(cnv_burden_chr_CNV_prop_broad)
  dim(cnv_burden_chr_CNV_prop_broad)
  ## [1] 46 26
  str(cnv_burden_chr_CNV_prop_broad)
  
  # 字符转数值
  mat <- as.data.frame(apply(cnv_burden_chr_CNV_prop_broad[,-1], 2, as.numeric))
  rownames(mat) <- cnv_burden_chr_CNV_prop_broad[,"HXID"]
  mat[1:10,]
  str(mat)
  dim(mat)
  ## [1] 46 25 44+2淋巴转移样本
  
  # 按照groups$Tumor筛选样本
  groups=read.csv("12_gistic/12.gistic.sample.group.csv",stringsAsFactors = F)
  groups$group[groups$group == "adenocarcinoma"] <- "PCA"
  groups$full_hxid <- paste(groups$hxid, groups$group, sep = "_")
  
  pos <- match(groups$full_hxid, rownames(mat))
  length(pos)## 44
  mat <- as.matrix(mat[pos, ])
  
  str(mat)
  mat[1:10,]
  dim(mat)
  ## [1] 44 22
 
  set.seed(13)
  col_ha =  columnAnnotation(Cancer = groups[,c('group')],
                             ADT = groups[,c("ADT")],
                             Mutation_group = groups[, c('Mutation_group')],
                             col=list(Mutation_group = c("A" = "#48D1CC", "B" = "#FFD700"),
                                      ADT = c("0" = brewer.pal(7,"Set1")[3], "1" = brewer.pal(7,"Set1")[4]),
                                      Cancer = c("IDCP" = brewer.pal(7,"Set1")[1] , "PCA" = brewer.pal(7,"Set1")[2]) ))
  
  f1 = colorRamp2(seq(min(mat), max(mat), length = 2), c("white", "red"))
  pdf(file.path(outdir, "heatmap_CNV_burden.pdf"), width = 8, height = 5)
  Heatmap(t(mat), # 注意转置
          name = "CNV burden",
          col = f1,
          cluster_rows = T,
          cluster_columns = T,
          bottom_annotation = col_ha,
          row_names_gp = gpar(fontsize = 11),
          column_names_gp = gpar(fontsize = 7)
  )
  dev.off()
  
  # 方差分析
  chr_aov <- aov(CNV_length_prop ~ Chr, cnv_burden_chr)
  summary(chr_aov)
  
  # 各水平均值
  print(model.tables(chr_aov,"means"),digits=3)
  
  # 多重t检验方法针对每一水平数据进行t检验。
  TukeyHSD(chr_aov)
  
  # https://www.liujason.com/article/155.html
  result_O <- HSD.test(chr_aov, "Chr", group = T)
  print(result_O)
  
  chr_aov
  tukey <- glht(chr_aov, linfct = mcp(Chr = "Tukey"))
  tukey
  #summary(tukey)
  tukey.cld <- cld(tukey)
  #opar <- par(mai=c(1,1,1.6,1))
  plot(tukey.cld)
  #par(opar)
}



#### 2. CNV burden: CNV gain length / chromosome length 
{
  # 使用spread函数将gd1_long长数据转换为宽数据gd1_wide
  cnv_burden_chr_gain_length_prop_broad <- spread(cnv_burden_chr[, c("HXID", "Chr", "Gian_length_prop")], Chr, Gian_length_prop)
  head(cnv_burden_chr_gain_length_prop_broad)
  dim(cnv_burden_chr_gain_length_prop_broad)
  ## [1] 46 23
  str(cnv_burden_chr_gain_length_prop_broad)
  
  # 字符转数值
  mat <- as.data.frame(apply(cnv_burden_chr_gain_length_prop_broad[,-1], 2, as.numeric))
  rownames(mat) <- cnv_burden_chr_gain_length_prop_broad[,"HXID"]
  mat[1:10,]
  dim(mat)
  ## [1] 46 22 44+2淋巴转移样本
  
  pos <- match(groups$full_hxid, rownames(mat))
  length(pos) ## 44
  mat <- as.matrix(mat[pos, ])
  dim(mat)
  
  set.seed(13)
  col_ha =  columnAnnotation(Cancer = groups[,c('group')],
                             ADT = groups[,c("ADT")],
                             Mutation_group = groups[, c('Mutation_group')],
                             col=list(Mutation_group = c("A" = "#48D1CC", "B" = "#FFD700"),
                                      ADT = c("0" = brewer.pal(7,"Set1")[3], "1" = brewer.pal(7,"Set1")[4]),
                                      Cancer = c("IDCP" = brewer.pal(7,"Set1")[1] , "PCA" = brewer.pal(7,"Set1")[2]) ))
  
  f1 = colorRamp2(seq(min(mat), max(mat), length = 2), c("white", "red"))
  pdf(file.path(outdir, "heatmap_CNV_burden_gain_length_proportion.pdf"), width = 8, height = 5)
  Heatmap(t(mat), # 注意转置
          name = "CNV burden",
          col = f1,
          cluster_rows = T,
          cluster_columns = T,
          bottom_annotation = col_ha,
          row_names_gp = gpar(fontsize = 11),
          column_names_gp = gpar(fontsize = 7)
  )
  dev.off()
  
  # 方差分析
  chr_aov <- aov(cnv_burden_chr$Gian_length_prop ~ Chr, cnv_burden_chr)
  summary(chr_aov)
  # 各水平均值
  print(model.tables(chr_aov,"means"),digits=3)
  # 多重t检验方法针对每一水平数据进行t检验。
  hsd <- TukeyHSD(chr_aov)
  hsd
#  plot(chr_aov)
  
  # https://www.liujason.com/article/155.html
  library(multcomp)
  tukey <- glht(chr_aov, linfct = mcp(Chr = "Tukey"))
  tukey
  #summary(tukey)
  tukey.cld <- cld(tukey)
  tukey.cld
  opar <- par(mai=c(1,1,1.6,1))
  plot(tukey.cld)
  par(opar)
}



#### 3. CNV burden: CNV loss length / chromosome length 
{
  # 使用spread函数将gd1_long长数据转换为宽数据gd1_wide
  cnv_burden_chr_loss_length_prop_broad <- spread(cnv_burden_chr[, c("HXID", "Chr", "Loss_length_prop")], Chr, Loss_length_prop)
  head(cnv_burden_chr_loss_length_prop_broad)
  dim(cnv_burden_chr_loss_length_prop_broad)
  ## [1] 46 23
  str(cnv_burden_chr_loss_length_prop_broad)
  
  # 字符转数值
  mat <- as.data.frame(apply(cnv_burden_chr_loss_length_prop_broad[,-1], 2, as.numeric))
  rownames(mat) <- cnv_burden_chr_loss_length_prop_broad[,"HXID"]
  mat[1:10,]
  dim(mat)
  ## [1] 46 22 44+2淋巴转移样本
  
  pos <- match(groups$full_hxid, rownames(mat))
  length(pos) ## 44
  mat <- as.matrix(mat[pos, ])
  dim(mat)
  
  set.seed(13)
  col_ha =  columnAnnotation(Cancer = groups[,c('group')],
                             ADT = groups[,c("ADT")],
                             Mutation_group = groups[, c('Mutation_group')],
                             col=list(Mutation_group = c("A" = "#48D1CC", "B" = "#FFD700"),
                                      ADT = c("0" = brewer.pal(7,"Set1")[3], "1" = brewer.pal(7,"Set1")[4]),
                                      Cancer = c("IDCP" = brewer.pal(7,"Set1")[1] , "PCA" = brewer.pal(7,"Set1")[2]) ))
  
  f1 = colorRamp2(seq(min(mat), max(mat), length = 2), c("white", "red"))
  pdf(file.path(outdir, "heatmap_CNV_burden_loss_length_proportion.pdf"), width = 8, height = 5)
  Heatmap(t(mat), # 注意转置
          name = "CNV burden",
          col = f1,
          cluster_rows = T,
          cluster_columns = T,
          bottom_annotation = col_ha,
          row_names_gp = gpar(fontsize = 11),
          column_names_gp = gpar(fontsize = 7)
  )
  dev.off()
  
  # 方差分析
  chr_aov <- aov(cnv_burden_chr$Loss_length_prop ~ Chr, cnv_burden_chr)
  summary(chr_aov)
  # 各水平均值
  print(model.tables(chr_aov,"means"),digits=3)
  # 多重t检验方法针对每一水平数据进行t检验。
  TukeyHSD(chr_aov)
  
  # https://www.liujason.com/article/155.html
  tukey <- glht(chr_aov, linfct = mcp(Chr = "Tukey"))
  tukey
  #summary(tukey)
  tukey.cld <- cld(tukey)
  tukey.cld
  opar <- par(mai=c(1,1,1.6,1))
  plot(tukey.cld)
  par(opar)
}


# 检查异常显著的
mat["43674S01X3",]
gain <- dat[which(dat$tumor == "43674S01X3" & dat$CNV_classification == "gain"), c("chromosome", "start.pos", "end.pos")]
sum(gain$end.pos - gain$start.pos + 1) / 63025520
## 59,590,497
dat[which(dat$tumor == "43674S01X3" & dat$CNV_classification == "loss"), c("chromosome", "start.pos", "end.pos")]






############################### Part 2 ############################### 
# 不区分染色体

cnv_burden=c()
for (sample in samples) {
  # 1) the total length of genomic DNA involved in identified CNV events
  sub_gain = dat[which(dat$full_hxid == sample & dat$CNV_classification == "gain"),]
  sub_gain_length = sum(sub_gain$end.pos - sub_gain$start.pos + 1)
  
  sub_loss = dat[which(dat$full_hxid == sample & dat$CNV_classification == "loss"),]
  sub_loss_length = sum(sub_loss$end.pos - sub_loss$start.pos + 1)
  
  sub_total_length <- sub_gain_length + sub_loss_length
  
  # 2) total number of CNVs
  sub = dat[dat$full_hxid == sample, ]
  sub_gain_number = nrow(sub[which(sub$CNV_classification == "gain"), ])
  sub_loss_number = nrow(sub[which(sub$CNV_classification == "loss"), ])
  sub_total_number = sub_gain_number + sub_loss_number
  
  # 2) average length of CNVs
  sub_gain_ave_len <- sub_gain_length / sub_gain_number
  sub_loss_ave_len <- sub_loss_length / sub_loss_number
  sub_total_ave_len  <- sub_total_length / sub_total_number
  
  cnv_burden=rbind(cnv_burden, c(sample, sub_gain_length, sub_loss_length, sub_total_length, 
                                 sub_gain_number, sub_loss_number, sub_total_number,
                                 sub_gain_ave_len, sub_loss_ave_len, sub_total_ave_len))
  
}
cnv_burden=as.data.frame(cnv_burden, stringsAsFactors =F)
colnames(cnv_burden)=c("HXID","Gain_length","Loss_length","Total_length", 
                       "Gain_number", "Loss_number", "Total_number", 
                       "Ave_length_gain", "Ave_length_loss", "Ave_length_total")
str(cnv_burden)
cnv_burden$Gain_length = as.numeric(cnv_burden$Gain_length)
cnv_burden$Loss_length = as.numeric(cnv_burden$Loss_length)
cnv_burden$Total_length = as.numeric(cnv_burden$Total_length)
cnv_burden$Gain_number <- as.numeric(cnv_burden$Gain_number)
cnv_burden$Loss_number <- as.numeric(cnv_burden$Loss_number)
cnv_burden$Total_number <- as.numeric(cnv_burden$Total_number)
cnv_burden$Ave_length_gain <- as.numeric(cnv_burden$Ave_length_gain)
cnv_burden$Ave_length_loss <- as.numeric(cnv_burden$Ave_length_loss)
cnv_burden$Ave_length_total <- as.numeric(cnv_burden$Ave_length_total)
cnv_burden[1:2, ]

groups

pos = match(groups$full_hxid, cnv_burden$HXID)
res = cbind(groups[,-6], cnv_burden[pos,])
#res$group[res$group=="adenocarcinoma"]="PCA"
head(res)
write.csv(res, file.path(outdir, "cnv_burden_20211214.csv"), quote = F, row.names = F)
dim(res)
## [1] 44 18


library(ggpubr)
####### IDCP vs PCA
my_comparisons <- list( c("IDCP", "PCA"))
p1=ggpaired(res, x = "group", y = "Gain_number", shape = "group", color = "group", id= "hxid", palette = "npg", line.color = "gray", xlab="", linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = T)+ylab("CNV gain number")
p1
p2=ggpaired(res, x = "group", y = "Loss_number",shape = "group",color = "group",id="hxid",palette = "npg",line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("CNV loss number")
p2
p3=ggpaired(res, x = "group", y = "Total_number",shape = "group",color = "group",id="hxid",palette = "npg",line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Total CNV number")
p3
ggarrange(p1+ylim(0,320),p2+ylim(0,60),p3+ylim(0,350),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.01.sequenza.cnv.burden.number.idcp_vs_pca.pdf"),width = 8,height = 4)

t.test(res$Gain_number ~ res$group, res)
t.test(res$Loss_number ~ res$group, res)
t.test(res$Total_number ~ res$group, res)


my_comparisons <- list( c("IDCP", "PCA"))
p1=ggpaired(res, x = "group", y = "Gain_length",shape = "group",color = "group",id="hxid",palette = "npg",line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("CNV gain length (bp)")
p1
p2=ggpaired(res, x = "group", y = "Loss_length",shape = "group",color = "group",id="hxid",palette = "npg",line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("CNV loss length (bp)")
p2
p3=ggpaired(res, x = "group", y = "Total_length",shape = "group",color = "group",id="hxid",palette = "npg",line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Total CNV length (bp)")
p3
ggarrange(p1+ylim(0, 2.3e9),p2+ylim(0, 1.35e9),p3+ylim(0, 2.8e9),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.02.sequenza.cnv.burden.length.idcp_vs_pca.pdf"),width = 8,height = 4)


my_comparisons <- list( c("IDCP", "PCA"))
p1=ggpaired(res, x = "group", y = "Ave_length_gain", shape = "group",color = "group",id="hxid",palette = "npg",line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Average length (bp) of CNV gain")
p1
p2=ggpaired(res, x = "group", y = "Ave_length_loss", shape = "group",color = "group",id="hxid",palette = "npg",line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Average length (bp) of CNV loss")
p2
p3=ggpaired(res, x = "group", y = "Ave_length_total", shape = "group",color = "group",id="hxid",palette = "npg",line.color = "gray",xlab="",linetype = "twodash")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("Average length (bp) of CNVs")
p3
ggarrange(p1+ylim(0,4e7),p2+ylim(0,4.5e7),p3+ylim(0,3e7), ncol=3, nrow=1, common.legend = T, labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.03.sequenza.cnv.burden.ave_length.idcp_vs_pca.pdf"),width = 8,height = 4)




####### ADT vs non_ADT
## CNV number
my_comparisons <- list( c("1", "0"))
res$ADT <- as.factor(res$ADT)
str(res)
p1=ggboxplot(res, x = "ADT", y = "Gain_number",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test", paired = F)+ylab("CNV gain number")
p1
p2=ggboxplot(res, x = "ADT", y = "Loss_number",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss number")
p2
p3=ggboxplot(res, x = "ADT", y = "Total_number",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV number")
p3
ggarrange(p1+ylim(0, 320),p2+ylim(0, 60),p3+ylim(0, 350),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.04.sequenza.cnv.burden.number.ADT_vs_nonADT.pdf"),width = 8,height = 4)

## CNV length
my_comparisons <- list( c("1", "0"))
p1=ggboxplot(res, x = "ADT", y = "Gain_length",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain length (bp)")
p1
p2=ggboxplot(res, x = "ADT", y = "Loss_length",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss length (bp)")
p2
p3=ggboxplot(res, x = "ADT", y = "Total_length",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV length (bp)")
p3
ggarrange(p1+ylim(0,2.3e9),p2+ylim(0,1.35e9),p3+ylim(0,2.8e9),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.05.sequenza.cnv.burden.length.ADT_vs_nonADT.pdf"),width = 8,height = 4)

# Average CNV length
my_comparisons <- list( c("1", "0"))
p1=ggboxplot(res, x = "ADT", y = "Ave_length_gain",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV gain")
p1
p2=ggboxplot(res, x = "ADT", y = "Ave_length_loss",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV loss")
p2
p3=ggboxplot(res, x = "ADT", y = "Ave_length_total",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNVs")
p3
ggarrange(p1+ylim(0, 4e7),p2+ylim(0, 4.5e7),p3+ylim(0,3e7),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.06.sequenza.cnv.burden.ave_length.ADT_vs_nonADT.pdf"),width = 8,height = 4)




####### A vs B
my_comparisons <- list( c("A", "B"))
p1=ggboxplot(res, x = "Mutation_group", y = "Gain_length",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain length (bp)")
p1
p2=ggboxplot(res, x = "Mutation_group", y = "Loss_length",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV oss length (bp)")
p2
p3=ggboxplot(res, x = "Mutation_group", y = "Total_length",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV length (bp)")
p3
ggarrange(p1+ylim(0,2.3e9),p2+ylim(0,1.35e9),p3+ylim(0,2.85e9),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.07.sequenza.cnv.burden.length.A_vs_B.pdf"),width = 8,height = 4)


p1=ggboxplot(res, x = "Mutation_group", y = "Gain_number",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain number")
p1
p2=ggboxplot(res, x = "Mutation_group", y = "Loss_number",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss number")
p2
p3=ggboxplot(res, x = "Mutation_group", y = "Total_number",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV number")
p3
ggarrange(p1+ylim(0,320), p2+ylim(0,60), p3+ylim(0,350), ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.08.sequenza.cnv.burden.number.A_vs_B.pdf"),width = 8,height = 4)



my_comparisons <- list( c("A", "B"))
p1=ggboxplot(res, x = "Mutation_group", y = "Ave_length_gain",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV gain")
p1
p2=ggboxplot(res, x = "Mutation_group", y = "Ave_length_loss",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV loss")
p2
p3=ggboxplot(res, x = "Mutation_group", y = "Ave_length_total",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNVs")
p3
ggarrange(p1+ylim(0,4e7),p2+ylim(0,4.5e7),p3+ylim(0,3e7),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.09.sequenza.cnv.burden.ave_length.A_vs_B.pdf"),width = 8,height = 4)




## within IDCP
idcp=res[res$group=="IDCP",]

####### within IDCP: ADT vs non_ADT
## number
my_comparisons <- list( c("1", "0"))

p1=ggboxplot(idcp, x = "ADT", y = "Gain_number",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain number within IDCP")
p1
p2=ggboxplot(idcp, x = "ADT", y = "Loss_number",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss number within IDCP")
p2
p3=ggboxplot(idcp, x = "ADT", y = "Total_number",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV number within IDCP")
p3
ggarrange(p1+ylim(0,320),p2+ylim(0,60),p3+ylim(0,350),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.10.sequenza.cnv.burden.number.within.idcp.ADT_vs_nonADT.pdf"),width = 8,height = 4)


####### within IDCP: ADT vs non_ADT
## length
p1=ggboxplot(idcp, x = "ADT", y = "Gain_length",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain length (bp) within IDCP")
p1
p2=ggboxplot(idcp, x = "ADT", y = "Loss_length",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss length (bp) within IDCP")
p2
p3=ggboxplot(idcp, x = "ADT", y = "Total_length",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV length (bp) within IDCP")
p3
ggarrange(p1+ylim(0,7.5e8),p2+ylim(0,1.45e9),p3+ylim(0,1.7e9),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.11.sequenza.cnv.burden.length.within.idcp.ADT_vs_nonADT.pdf"),width = 8,height = 4)

####### within IDCP: ADT vs non_ADT
####### ADT vs non ADT
## ave length
p1=ggboxplot(idcp, x = "ADT", y = "Ave_length_gain",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV gain within IDCP")
p1
p2=ggboxplot(idcp, x = "ADT", y = "Ave_length_loss",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV loss within IDCP")
p2
p3=ggboxplot(idcp, x = "ADT", y = "Ave_length_total",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNVs within IDCP")
p3
ggarrange(p1+ylim(0, 4e7),p2+ylim(0, 4.5e7),p3+ylim(0, 3e7),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.12.sequenza.cnv.burden.ave_length.within.idcp.ADT_vs_nonADT.pdf"),width = 8,height = 4)



####### within IDCP: ADT vs non_ADT
####### A vs B
my_comparisons <- list( c("A", "B"))

####### within IDCP: ADT vs non_ADT
# Number 
my_comparisons <- list( c("A", "B"))
p1=ggboxplot(idcp, x = "Mutation_group", y = "Gain_number",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain number within IDCP")
p1
p2=ggboxplot(idcp, x = "Mutation_group", y = "Loss_number",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss number within IDCP")
p2
p3=ggboxplot(idcp, x = "Mutation_group", y = "Total_number",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV number within IDCP")
p3
ggarrange(p1+ylim(0,320),p2+ylim(0,60),p3+ylim(0,350),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.13.sequenza.cnv.burden.number.within.idcp.A_vs_B.pdf"),width = 8,height = 4)


####### within IDCP: ADT vs non_ADT
# length
my_comparisons <- list( c("A", "B"))
p1=ggboxplot(idcp, x = "Mutation_group", y = "Gain_length",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain length (bp) within IDCP")
p1
p2=ggboxplot(idcp, x = "Mutation_group", y = "Loss_length",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss length (bp) within IDCP")
p2
p3=ggboxplot(idcp, x = "Mutation_group", y = "Total_length",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV length (bp) within IDCP")
p3
ggarrange(p1+ylim(0,7.6e8),p2+ylim(0,1.4e9),p3+ylim(0,1.6e9),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.14.sequenza.cnv.burden.length.within.idcp.A_vs_B.pdf"),width = 8,height = 4)


####### within IDCP: ADT vs non_ADT
# average length
my_comparisons <- list( c("A", "B"))
p1=ggboxplot(idcp, x = "Mutation_group", y = "Ave_length_gain",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNA gain within IDCP")
p1
p2=ggboxplot(idcp, x = "Mutation_group", y = "Ave_length_loss",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV loss within IDCP")
p2
p3=ggboxplot(idcp, x = "Mutation_group", y = "Ave_length_total",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNVs within IDCP")
p3
ggarrange(p1+ylim(0, 4e7),p2+ylim(0, 4.5e7),p3+ylim(0, 3e7),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.15.sequenza.cnv.burden.ave_length.within.idcp.A_vs_B.pdf"),width = 8,height = 4)




## within pca
pca=res[res$group=="PCA",]
####### within pca: ADT vs non_ADT
## Number 
my_comparisons <- list( c("1", "0"))
p1=ggboxplot(pca, x = "ADT", y = "Gain_number",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain number within PCA")
p1
p2=ggboxplot(pca, x = "ADT", y = "Loss_number",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss number within PCA")
p2
p3=ggboxplot(pca, x = "ADT", y = "Total_number",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV number within PCA")
p3
ggarrange(p1+ylim(0, 240),p2+ylim(0, 55),p3+ylim(0, 260),ncol=3,nrow=1, common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.16.sequenza.cnv.burden.number.within.pca.ADT_vs_nonADT.pdf"),width = 8,height = 4)


####### within pca: ADT vs non_ADT
## length 
my_comparisons <- list( c("1", "0"))
p1=ggboxplot(pca, x = "ADT", y = "Gain_length",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain length (bp) within PCA")
p1
p2=ggboxplot(pca, x = "ADT", y = "Loss_length",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss length (bp) within PCA")
p2
p3=ggboxplot(pca, x = "ADT", y = "Total_length",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total length (bp) within PCA")
p3
ggarrange(p1+ylim(0,2.3e9),p2+ylim(0,7.5e8),p3+ylim(0,2.8e9),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.17.sequenza.cnv.burden.length.within.pca.ADT_vs_nonADT.pdf"),width = 8,height = 4)


####### within pca: ADT vs non_ADT
## ave length 
my_comparisons <- list( c("1", "0"))
p1=ggboxplot(pca, x = "ADT", y = "Ave_length_gain",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV gain within PCA")
p1
p2=ggboxplot(pca, x = "ADT", y = "Ave_length_loss",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV loss within PCA")
p2
p3=ggboxplot(pca, x = "ADT", y = "Ave_length_total",shape = "ADT",color = "ADT",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNVs within PCA")
p3
ggarrange(p1+ylim(0, 2.1e7),p2+ylim(0, 3.6e7),p3+ylim(0, 2e7), ncol=3, nrow=1, common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.18.sequenza.cnv.burden.ave_length.within.pca.ADT_vs_nonADT.pdf"),width = 8,height = 4)




####### within pca: A vs B
####### A vs B
## Number
my_comparisons <- list( c("A", "B"))
p1=ggboxplot(pca, x = "Mutation_group", y = "Gain_number",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain number within PCA")
p1
p2=ggboxplot(pca, x = "Mutation_group", y = "Loss_number",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss number within PCA")
p2
p3=ggboxplot(pca, x = "Mutation_group", y = "Total_number",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV number within PCA")
p3
ggarrange(p1+ylim(0, 240),p2+ylim(0, 55),p3+ylim(0, 260),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.19.sequenza.cnv.burden.number.within.pca.A_vs_B.pdf"),width = 8,height = 4)


####### within pca: A vs B
## length
my_comparisons <- list( c("A", "B"))
p1=ggboxplot(pca, x = "Mutation_group", y = "Gain_length",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV gain length (bp) within PCA")
p1
p2=ggboxplot(pca, x = "Mutation_group", y = "Loss_length",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("CNV loss length (bp) within PCA")
p2
p3=ggboxplot(pca, x = "Mutation_group", y = "Total_length",shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Total CNV length (bp) within PCA")
p3
ggarrange(p1+ylim(0,2.3e9),p2+ylim(0,7.5e8),p3+ylim(0,2.8e9),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.20.sequenza.cnv.burden.length.within.pca.A_vs_B.pdf"),width = 8,height = 4)

pca[1:2,]
pca
t.test(pca$Gain_length~pca$Mutation_group, pca)
t.test(pca$Loss_length~pca$Mutation_group, pca)
t.test(pca$Total_length~pca$Mutation_group, pca)


####### within pca: A vs B
## ave length 
my_comparisons <- list( c("A", "B"))
p1=ggboxplot(pca, x = "Mutation_group", y = "Ave_length_gain", shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV gain within PCA")
p1
p2=ggboxplot(pca, x = "Mutation_group", y = "Ave_length_loss", shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNV loss within PCA")
p2
p3=ggboxplot(pca, x = "Mutation_group", y = "Ave_length_total" ,shape = "Mutation_group",color = "Mutation_group",palette = "npg",xlab="",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = F)+ylab("Average length (bp) of CNVs within PCA")
p3
ggarrange(p1+ylim(0, 2.1e7),p2+ylim(0,3.5e7),p3+ylim(0,1.9e7),ncol=3,nrow=1,common.legend = T,labels = c("A)","B)","C)"))
ggsave(file.path(outdir,"07.21.sequenza.cnv.burden.ave_length.within.pca.A_vs_B.pdf"),width = 8,height = 4)




