
# 1. 使用maftools对goistic的结果进行可视化
# 2. 比较IDCP和PCA之间差异的Gscore，gistic的G score均是正值，通过amp和del区分方向，但是maftools里面对G score直接区分了正负
#    由于亚群间的segement是不同的，取决于各自亚群的特点，所以无法给予segment直接进行比较，只能对两个亚群间存在oevrlap segment进行G socre 的比较
#    同时将amp/del转换为+/-，转移到G scores上，这样可以直接计算diff G socre
#    IDCP_vs_PCA diff_G_score > 0: IDCP amp 大于PCA  

# 2021-12-10: 
## 检查IDCP vs PCA cytoband Sample frequency heat-map, 发现了unique错了导致的heatmap里面缺失出现问题，结果


rm(list = ls())
setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ/")

library(dplyr)
library(ComplexHeatmap)
library(maftools)

# 根据运行结果，依次将结果赋值给相应的参数，以构建GISTIC输入文件。
all44_chrxy.gistic <- readGistic(gisticAllLesionsFile = '12.sequenza.44.samples.gistic.chrx/all_lesions.conf_90.txt',
                          gisticAmpGenesFile = '12.sequenza.44.samples.gistic.chrx/amp_genes.conf_90.txt',
                          gisticDelGenesFile = '12.sequenza.44.samples.gistic.chrx/del_genes.conf_90.txt',
                          gisticScoresFile = '12.sequenza.44.samples.gistic.chrx/scores.gistic',
                          isTCGA = F)

# 接下来，对结果进行可视化展示
pdf("12.sequenza.44.samples.gistic.chrx/all44_chrxy_gisticChromPlot.pdf", width = 7, height = 3)
gisticChromPlot(gistic = all44_chrxy.gistic, ref.build = 'hg19', fdrCutOff = 0.1, cytobandTxtSize = 0.8)
# 结果显示：其中红色表示拷贝数增加，而蓝色表示拷贝数降低，并对其中几个较为显著的基因名字信息进行了标注。
dev.off()

pdf("12.sequenza.44.samples.gistic.chrx/all44_chrxy_gisticBubblePlot.pdf", width = 5.5, height = 5)
gisticBubblePlot(gistic = all44_chrxy.gistic, fdrCutOff = 0.1)
dev.off()

#oncoplot
pdf(file.path("12.sequenza.44.samples.gistic.chrx/", "all44_chrxy_gisticOncoPlot.pdf"), height = 12, width = 7)
gisticOncoPlot(gistic = all44_chrxy.gistic, 
               showTumorSampleBarcodes = T,
               sepwd_genes = 1,
               sepwd_samples = 1,
               SampleNamefontSize = 0.7,
               fontSize = 0.7,
               legendFontSize = 1,
               borderCol = "black",
               bgCol = "white",
               colors = c("Amp" = "red", "Del" = "blue"),
               sortByAnnotation = T)
dev.off()






################################################################################################################################
# 分开IDCP和PCA样本运行gistic2

# IDCP 
rm(list = ls())
setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ/")

library(maftools)
# 根据运行结果，依次将结果赋值给相应的参数，以构建GISTIC输入文件。
# Linux 
IDCP_chrxy.gistic <- readGistic(gisticAllLesionsFile = '12.sequenza.44.samples.gistic.chrx.IDCP/all_lesions.conf_90.txt',
                                 gisticAmpGenesFile = '12.sequenza.44.samples.gistic.chrx.IDCP/amp_genes.conf_90.txt',
                                 gisticDelGenesFile = '12.sequenza.44.samples.gistic.chrx.IDCP/del_genes.conf_90.txt',
                                 gisticScoresFile = '12.sequenza.44.samples.gistic.chrx.IDCP/scores.gistic',
                                 isTCGA = F)

# GenePatern的GISTIC2的结果
#IDCP_chrxy.gistic <- readGistic(gisticAllLesionsFile = '12.sequenza.44.samples.gistic.chrx.IDCP/GenePatern_Gistic/394095/all_lesions.conf_90.txt',
#                                gisticAmpGenesFile = '12.sequenza.44.samples.gistic.chrx.IDCP/GenePatern_Gistic/394095/amp_genes.conf_90.txt',
#                                gisticDelGenesFile = '12.sequenza.44.samples.gistic.chrx.IDCP/GenePatern_Gistic/394095/del_genes.conf_90.txt',
#                                gisticScoresFile = '12.sequenza.44.samples.gistic.chrx.IDCP/GenePatern_Gistic/394095/scores.gistic',
#                                isTCGA = F)



# Linux
PCA_chrxy.gistic <- readGistic(gisticAllLesionsFile = '12.sequenza.44.samples.gistic.chrx.PCA/all_lesions.conf_90.txt',
                               gisticAmpGenesFile = '12.sequenza.44.samples.gistic.chrx.PCA/amp_genes.conf_90.txt',
                               gisticDelGenesFile = '12.sequenza.44.samples.gistic.chrx.PCA/del_genes.conf_90.txt',
                               gisticScoresFile = '12.sequenza.44.samples.gistic.chrx.PCA/scores.gistic',
                               isTCGA = F)

# GenePatern的GISTIC2的结果
#PCA_chrxy.gistic <- readGistic(gisticAllLesionsFile = '12.sequenza.44.samples.gistic.chrx.PCA/GenePattern_Gistic/394096/all_lesions.conf_90.txt',
#                               gisticAmpGenesFile = '12.sequenza.44.samples.gistic.chrx.PCA/GenePattern_Gistic/394096/amp_genes.conf_90.txt',
#                               gisticDelGenesFile = '12.sequenza.44.samples.gistic.chrx.PCA/GenePattern_Gistic/394096/del_genes.conf_90.txt',
#                               gisticScoresFile = '12.sequenza.44.samples.gistic.chrx.PCA/GenePattern_Gistic/394096/scores.gistic',
#                               isTCGA = F)


PCA_chrxy.gistic@data
PCA_chrxy.gistic@gis.scores


# 接下来，对结果进行可视化展示
# IDCP
pdf(file.path("12.sequenza.44.samples.gistic.chrx.IDCP/", "IDCP_chrxy_gisticChromPlot.pdf"), width = 7, height = 3)
gisticChromPlot(gistic = IDCP_chrxy.gistic, ref.build = 'hg19', fdrCutOff = 0.1, cytobandTxtSize = 0.8)
# 结果显示：其中红色表示拷贝数增加，而蓝色表示拷贝数降低，并对其中几个较为显著的基因名字信息进行了标注。
dev.off()


pdf(file.path("12.sequenza.44.samples.gistic.chrx.IDCP/", "IDCP_chrxy_gisticBubblePlot.pdf"), width = 5.5, height = 5)
gisticBubblePlot(gistic = IDCP_chrxy.gistic)
dev.off()


#oncoplot
pdf(file.path("12.sequenza.44.samples.gistic.chrx.IDCP/", "IDCP_chrxy_gisticOncoPlot.pdf"), height = 9, width = 5)
gisticOncoPlot(gistic = IDCP_chrxy.gistic, 
               showTumorSampleBarcodes = T,
               sepwd_genes = 1,
               sepwd_samples = 1,
               SampleNamefontSize = 0.7,
               fontSize = 0.7,
               legendFontSize = 1,
               borderCol = "black",
               bgCol = "white",
               colors = c("Amp" = "red", "Del" = "blue"),
               sortByAnnotation = F)
dev.off()




############################################################################################################
# PCA

pdf(file.path("12.sequenza.44.samples.gistic.chrx.PCA/", "PCA_chrxy_gisticChromPlot.pdf"), width = 7, height = 3)
gisticChromPlot(gistic = PCA_chrxy.gistic, ref.build = 'hg19', fdrCutOff = 0.1, cytobandTxtSize = 0.8)
# 结果显示：其中红色表示拷贝数增加，而蓝色表示拷贝数降低，并对其中几个较为显著的基因名字信息进行了标注。
dev.off()

pdf(file.path("12.sequenza.44.samples.gistic.chrx.PCA/", "PCA_chrxy_gisticBubblePlot.pdf"), width = 5.5, height = 5)
gisticBubblePlot(gistic = PCA_chrxy.gistic)
dev.off()

#oncoplot
pdf(file.path("12.sequenza.44.samples.gistic.chrx.PCA/", "PCA_chrxy_gisticOncoPlot.pdf"), height = 9, width = 5)
gisticOncoPlot(gistic = PCA_chrxy.gistic, 
               showTumorSampleBarcodes = T,
               sepwd_genes = 1,
               sepwd_samples = 1,
               SampleNamefontSize = 0.7,
               fontSize = 0.7,
               legendFontSize = 1,
               borderCol = "black",
               bgCol = "white",
               colors = c("Amp" = "red", "Del" = "blue"),
               sortByAnnotation = T)
dev.off()

str(PCA_chrxy.gistic@cytoband.summary)









##############################################################################################################
# cytoband水平上比较IDCP和PCA
# q value < 0.25
PCA_cytoband <- as.data.frame(PCA_chrxy.gistic@cytoband.summary)
dim(PCA_cytoband)
## [1] 57  7
summary(PCA_cytoband)
write.csv(as.data.frame(PCA_chrxy.gistic@cytoband.summary), 
          file.path("12.sequenza.44.samples.gistic.chrx.PCA/", "Heatmap_PCA_cytoband_sunmary_info.csv"), quote = F, row.names = F)


IDCP_cytoband <- as.data.frame(IDCP_chrxy.gistic@cytoband.summary)
dim(IDCP_cytoband)
## [1] 52  7
summary(IDCP_cytoband)
head(IDCP_cytoband,10)
write.csv(as.data.frame(IDCP_chrxy.gistic@cytoband.summary), 
          file.path("12.sequenza.44.samples.gistic.chrx.IDCP/", "Heatmap_IDCP_cytoband_sunmary_info.csv"), quote = F, row.names = F)


# 合并cytoband后进行排序比较
a <- PCA_cytoband
a$Cancer <- "PCA"

b <- IDCP_cytoband
b$Cancer <- "IDCP"

IDCP_PCA_cytoband <- rbind(a, b)
head(IDCP_PCA_cytoband)
dim(IDCP_PCA_cytoband)
## [1] 109   8 57+52 = 109 
dim(unique(IDCP_PCA_cytoband))
## [1] 109   8
IDCP_PCA_cytoband <- IDCP_PCA_cytoband[order(IDCP_PCA_cytoband$Cytoband, decreasing = T),]
IDCP_PCA_cytoband[which(IDCP_PCA_cytoband$Unique_Name == "13q14.2"),]


IDCP_cytoband$Unique_Name <- lapply(IDCP_cytoband$Unique_Name, function(x) strsplit(x, ":")[[1]][2])
IDCP_cytoband[1:10,]
dim(IDCP_cytoband)
## 52 7 
PCA_cytoband$Unique_Name <- lapply(PCA_cytoband$Unique_Name, function(x) strsplit(x, ":")[[1]][2])
PCA_cytoband[1:10,]
dim(PCA_cytoband)
## 57 7 

FullJoin_IDCP_PCA_cytoband <- IDCP_cytoband  %>% full_join(PCA_cytoband, by = "Unique_Name") 
dim(FullJoin_IDCP_PCA_cytoband)
## [1] 90  13
FullJoin_IDCP_PCA_cytoband[1:3,]

length(unique(FullJoin_IDCP_PCA_cytoband$Unique_Name))

FullJoin_IDCP_PCA_cytoband$Unique_Name <- as.factor(FullJoin_IDCP_PCA_cytoband$Unique_Name)

FullJoin_IDCP_PCA_cytoband_SampleFrq <- FullJoin_IDCP_PCA_cytoband %>% mutate(IDCP_freq = if_else(Variant_Classification.x == "Amp", nSamples.x/22, -(nSamples.x/22))) %>% 
  mutate(PCA_freq = if_else(Variant_Classification.y == "Amp", nSamples.y/22, -(nSamples.y/22))) %>% mutate_at(c("IDCP_freq", "PCA_freq"), ~replace(., is.na(.), 0) )
FullJoin_IDCP_PCA_cytoband_SampleFrq[1:15,]
FullJoin_IDCP_PCA_cytoband_SampleFrq$Unique_Name <- as.character(FullJoin_IDCP_PCA_cytoband_SampleFrq$Unique_Name)
colnames(FullJoin_IDCP_PCA_cytoband_SampleFrq)
## [1] "Unique_Name"              "nGenes.x"                 "nSamples.x"               "Variant_Classification.x"
## [5] "Cytoband.x"               "Wide_Peak_Limits.x"       "qvalues.x"                "nGenes.y"                
## [9] "nSamples.y"               "Variant_Classification.y" "Cytoband.y"               "Wide_Peak_Limits.y"      
## [13] "qvalues.y"                "IDCP_freq"                "PCA_freq"      

str(FullJoin_IDCP_PCA_cytoband_SampleFrq)
FullJoin_IDCP_PCA_cytoband_SampleFrq[which(FullJoin_IDCP_PCA_cytoband_SampleFrq$Unique_Name == "13q14.11"),]
FullJoin_IDCP_PCA_cytoband_SampleFrq[which(FullJoin_IDCP_PCA_cytoband_SampleFrq$Unique_Name == "4p11"),]


dim(FullJoin_IDCP_PCA_cytoband_SampleFrq)
## [1] 90 15

# 去除重复
FullJoin_IDCP_PCA_cytoband_SampleFrq <- unique(FullJoin_IDCP_PCA_cytoband_SampleFrq)
dim(FullJoin_IDCP_PCA_cytoband_SampleFrq)
## [1] 90 15

FullJoin_IDCP_PCA_cytoband_SampleFrq_uniq <- unique(FullJoin_IDCP_PCA_cytoband_SampleFrq[, c(1, 14, 15)])
head(FullJoin_IDCP_PCA_cytoband_SampleFrq_uniq)
dim(FullJoin_IDCP_PCA_cytoband_SampleFrq_uniq)
## [1] 85  3
FullJoin_IDCP_PCA_cytoband_SampleFrq_uniq[which(FullJoin_IDCP_PCA_cytoband_SampleFrq_uniq$Unique_Name == "13q14.11"),]


mat <- as.matrix(FullJoin_IDCP_PCA_cytoband_SampleFrq_uniq[,c(2, 3)])
rownames(mat) <- FullJoin_IDCP_PCA_cytoband_SampleFrq_uniq[,1]
summary(mat)
dim(mat) 
# [1] 85  2
mat[1:5,]

set.seed(123)
ComplexHeatmap::Heatmap(mat, name = "% Samples",border = F, show_column_dend = F,
                        rect_gp = gpar(col = "black", lwd =1),
                        clustering_distance_rows = function(m) dist(m),
                        #clustering_distance_rows = function(x, y) 1 - cor(x, y),
                        row_names_side = "right",
                        row_names_gp = gpar(fontsize = 8, col = "black"),
                        column_names_gp = gpar(fontsize = 10, col = "black"),
                        km = 6
                        )

set.seed(123)
pdf(file.path("12.sequenza.44.samples.gistic.chrx/", "Heatmap_cytoband_IDCP_PCA_SampleFrequency.pdf"), width = 3, height = 8.5)
ComplexHeatmap::Heatmap(mat, name = "% Samples",border = F, show_column_dend = F,
                        rect_gp = gpar(col = "black", lwd =1),
                        clustering_distance_rows = function(m) dist(m),
                        #clustering_distance_rows = function(x, y) 1 - cor(x, y),
                        row_names_side = "right",
                        row_names_gp = gpar(fontsize = 7, col = "black"),
                        column_names_gp = gpar(fontsize = 10, col = "black"),
                        km = 7
)
dev.off()











#########################################################################################################
# gene level

sum(IDCP_chrxy.gistic@cytoband.summary$nGenes)
## [1] 4491
sum(PCA_chrxy.gistic@cytoband.summary$nGenes)
## [1] 1233

PCA_gene <- as.data.frame(PCA_chrxy.gistic@gene.summary)
PCA_gene[1:10,]
dim(PCA_gene)
## [1] 1212    5
colnames(PCA_gene)[2:5] <- paste("PCA", colnames(PCA_gene[2:5]), sep = "_") 
PCA_gene$Cancer <- "PCA"


IDCP_gene <- as.data.frame(IDCP_chrxy.gistic@gene.summary)
IDCP_gene[1:10,]
dim(IDCP_gene)
## [1] 4352    5
colnames(IDCP_gene)[2:5] <- paste("IDCP", colnames(IDCP_gene[2:5]), sep = "_") 
IDCP_gene$Cancer <- "IDCP"

FullJoin_IDCP_PCA_gene <- IDCP_gene  %>% full_join(PCA_gene, by = "Hugo_Symbol") 
dim(FullJoin_IDCP_PCA_gene)
## [1] 4982    11
head(FullJoin_IDCP_PCA_gene)

FullJoin_IDCP_PCA_gene$Cancer.y <- "PCA"

# NA替换0
FullJoin_IDCP_PCA_gene_NA_replaced_0 <- FullJoin_IDCP_PCA_gene %>% mutate(replace(., is.na(.), 0) ) 
dim(FullJoin_IDCP_PCA_gene_NA_replaced_0)
## [1] 4982   11
head(FullJoin_IDCP_PCA_gene_NA_replaced_0)


mat <- as.matrix(FullJoin_IDCP_PCA_gene_NA_replaced_0[,c("IDCP_Amp", "IDCP_Del", "PCA_Amp", "PCA_Del")])
rownames(mat) <- FullJoin_IDCP_PCA_gene_NA_replaced_0[,1]
mat[,3] <- -mat[,3] 
mat[,4] <- -mat[,4]
summary(mat)

# gene太多了
#ComplexHeatmap::Heatmap(mat, name = "Amp/Del", border = F, show_column_dend = F,
                        #rect_gp = gpar(col = "black", lwd =1),
                        clustering_distance_rows = function(m) dist(m),
                        #clustering_distance_rows = function(x, y) 1 - cor(x, y),
                        row_names_side = "right",
                        row_names_gp = gpar(fontsize = 8, col = "black"),
                        column_names_gp = gpar(fontsize = 10, col = "black")
                        #km = 7
)


# 删除缺失值
na_rm_FullJoin_IDCP_PCA_gene <- na.omit(FullJoin_IDCP_PCA_gene)
dim(na_rm_FullJoin_IDCP_PCA_gene)
## [1] 582  11
head(na_rm_FullJoin_IDCP_PCA_gene)



mat <- na_rm_FullJoin_IDCP_PCA_gene[ ,c("IDCP_Amp", "PCA_Amp", "IDCP_Del", "PCA_Del")]
rownames(mat) <- na_rm_FullJoin_IDCP_PCA_gene[ ,1]
str(mat)
mat$IDCP_Del <- -mat$IDCP_Del
mat$PCA_Del <- -mat$PCA_Del

# https://www.bioinfo-scrounger.com/archives/122/
library(circlize)
library(dendextend)
# color for dendrogram
dend = hclust(dist(mat))
dend = color_branches(dend, k = 7)

ComplexHeatmap::Heatmap(mat, name = "Amp/Del", border = F, 
                        show_column_dend = F, 
                        cluster_columns = F,
                        cluster_rows = dend,
                        #rect_gp = gpar(col = "black", lwd =1),
                        #clustering_distance_rows = function(m) dist(m),
                        #clustering_distance_rows = function(x, y) 1 - cor(x, y),
                        row_names_side = "right",
                       # col = colorRamp2(c(0,22), c("white","red")),
                        row_names_gp = gpar(fontsize = 3, col = "black"),
                        column_names_gp = gpar(fontsize = 10, col = "black"),
                        #km = 7
                        )



# 合并比较
pdf(file.path("12.sequenza.44.samples.gistic.chrx/", "Compared_IDCP_PCA_chrxy_gisticChromPlot.pdf"), width = 8, height = 5)
par(mfrow=c(2,1))
gisticChromPlot(gistic = IDCP_chrxy.gistic, ref.build = 'hg19', fdrCutOff = 0.1, cytobandTxtSize = 0.8)
gisticChromPlot(gistic = PCA_chrxy.gistic, ref.build = 'hg19', fdrCutOff = 0.1, cytobandTxtSize = 0.8)
dev.off()



#layout(matrix(1:4,2,2),widths=c(3,1),heights=c(1,1))  #将当时装置按照m进行划分，宽度之比为3:1，高度之比为1:1
#layout.show(4)
#dev.list()   # 列出打开装置的列表，
#gisticChromPlot(gistic = IDCP_chrxy.gistic, ref.build = 'hg19', fdrCutOff = 0.1, cytobandTxtSize = 0.8)
#gisticChromPlot(gistic = PCA_chrxy.gistic, ref.build = 'hg19', fdrCutOff = 0.1, cytobandTxtSize = 0.8)
#gisticBubblePlot(gistic = IDCP_chrxy.gistic)
#gisticBubblePlot(gistic = PCA_chrxy.gistic)








##################################################################################################################
# https://www.biostars.org/p/448652/
#  using GISTIC2 to identify recurrent SCNAs (somatic copy number alterations) on a cohort of primary tumour samples matched with metastatic tumour samples.
# ran GISTIC2 separately on the primary tumours and then metastatic tumours. 
# For each recurrent amplification and deletion identified, there is a G-score.
# The G-score is calculated based on the amplitude of the aberration as well as the frequency of its occurrence across samples.

# Now for each amplification and deletion identified in the metastatic tumours, 
# I would like to compare the G-score to the G-score in the primary tumours to identify regions distinct to mets.

head(IDCP_chrxy.gistic@gis.scores)

library(GenomicRanges)

IDCP <- IDCP_chrxy.gistic
IDCP@cnv.summary # 不做CNV summary,这个结果应该宋博之前做过了
IDCP@cytoband.summary
IDCP@gene.summary
IDCP@summary


PCA <- PCA_chrxy.gistic
PCA@cnv.summary
PCA@cytoband.summary
PCA@gene.summary
PCA@summary


colnames(IDCP@gis.scores)
colnames(IDCP@gis.scores) <- paste0("IDCP_", colnames(IDCP@gis.scores))

head(IDCP@gis.scores)
str(IDCP@gis.scores)

IDCP.gr <- makeGRangesFromDataFrame(
  df = IDCP@gis.scores,
  keep.extra.columns = TRUE,
  seqnames.field = "IDCP_Chromosome",
  start.field = "IDCP_Start_Position",
  end.field = "IDCP_End_Position"
)
length(IDCP.gr@ranges)
## [1] 7525
head(IDCP.gr)



PCA@gis.scores

colnames(PCA@gis.scores) <- paste0("PCA_", colnames(PCA@gis.scores))
PCA.gr <- makeGRangesFromDataFrame(
  df = PCA@gis.scores,
  keep.extra.columns = TRUE,
  seqnames.field = "PCA_Chromosome",
  start.field = "PCA_Start_Position",
  end.field = "PCA_End_Position"
)
str(PCA.gr)
length(PCA.gr@ranges)
## [1] 7111

# find IDCP regions also in PCA that overlap
r_idcp <- subsetByOverlaps (IDCP.gr, PCA.gr, type = "equal")
str(r_idcp)
head(r_idcp)
length(r_idcp@ranges)
## [1] 246

# find IDCP regions also in PCA that overlap # subsetByOverlaps，它能够直接输出我们想要得到的重叠的ranges：
ol.gr  <- subsetByOverlaps(IDCP.gr, PCA.gr, type = "equal")
ol.gr

# index for PCA hits in IDCP
both.ind <- findOverlaps(query = IDCP.gr, subject = PCA.gr, type = "equal")

# grab PCA meta cols by index 
tmp.gr  <- PCA.gr[subjectHits(both.ind),]
tmp.gr

# add IDCP cols to ol.gr
mcols(ol.gr)[c(6:10)] <- mcols(tmp.gr)
head(ol.gr,n=2)


# In terms of calculating the difference in G scores at overlapping regions, does taking the log2 fold change IDCP over PCA seem like a good strategy?
# calc log2 FC
ol.gr$diff_gscore <- ol.gr$IDCP_G_Score - ol.gr$PCA_G_Score
head(ol.gr)
length(ol.gr@ranges)
## [1] 246
table(ol.gr@seqnames)
str(ol.gr)

ol.gr.df <- as.data.frame(ol.gr)
ol.gr.df$seqnames <- paste("chr", ol.gr.df$seqnames, sep = "")
head(ol.gr.df)
dim(ol.gr.df)
## [1] 246  16


# 两个overlap区间内的gain loss相反
dim(ol.gr.df[which(ol.gr.df$IDCP_Variant_Classification != ol.gr.df$PCA_Variant_Classification),])
# [1]  2 16
ol.gr.df[which(ol.gr.df$IDCP_Variant_Classification != ol.gr.df$PCA_Variant_Classification),]


# 检查Del的amplitude理论上不应该超过2？
# Amp的amplitude理论上不低于0？
tapply(ol.gr.df$IDCP_Avg_amplitude, INDEX=ol.gr.df$IDCP_Variant_Classification, FUN=summary)
tapply(ol.gr.df$PCA_Avg_amplitude, INDEX=ol.gr.df$PCA_Variant_Classification, FUN=summary)


# 筛选统计显著的结果， 0.25是gistic的的默认阈值
ol.gr.df.FDR <- ol.gr.df  %>% filter(IDCP_fdr < 0.25 & PCA_fdr < 0.25)
dim(ol.gr.df.FDR)
## [1] 160  16

# ol.gr.df.FDR <- ol.gr.df %>% filter(IDCP_fdr < 0.25 | PCA_fdr < 0.25)
# dim(ol.gr.df.FDR)
## [1] 178  16


# Gain和loss相反的区间
dim(ol.gr.df.FDR[which(ol.gr.df.FDR$IDCP_Variant_Classification != ol.gr.df.FDR$PCA_Variant_Classification), ])
# [1]  2 16
ol.gr.df.FDR[which(ol.gr.df.FDR$IDCP_Variant_Classification != ol.gr.df.FDR$PCA_Variant_Classification), ]

# 热图展示IDCP VS PCA的overlap G score
ol.gr.df.FDR[1:5, ]


tapply(ol.gr.df.FDR$diff_gscore, INDEX = ol.gr.df.FDR$IDCP_Variant_Classification, summary)


# segment overlapped region length
sum(ol.gr.df.FDR$width)/1000000
## [1] 32.32591 Mb

df <- ol.gr.df.FDR %>% mutate(Pos= paste(seqnames,start,end, sep="_")) %>% mutate(IDCP_G = ifelse(IDCP_Variant_Classification == "Amp", IDCP_G_Score, -IDCP_G_Score)) %>%
  mutate(PCA_G = ifelse(PCA_Variant_Classification == "Amp", PCA_G_Score, -PCA_G_Score)) %>% mutate(Diff_G = IDCP_G - PCA_G)
df <- df[order(df$Pos),]
head(df, 10)
dim(df)
## [1] 160  20

write.csv(df, file.path("12.sequenza.44.samples.gistic.chrx/", "Gscore_overlapped_segment_regions_IDCP_PCA.csv"), quote = F, row.names = F)



mat <- as.matrix(df[ ,c("IDCP_G", "PCA_G", "Diff_G")])
rownames(mat) <- df[,"Pos"]
summary(mat)
dim(mat) # [1] 160  3
mat[1:10,]
rownames(mat)
summary(mat)

set.seed(123)
ComplexHeatmap::Heatmap(mat, name = "G score",border = F, show_column_dend = F,
                        rect_gp = gpar(col = "black", lwd =1),
                        #clustering_distance_rows = function(m) dist(m),
                        #clustering_distance_rows = function(x, y) 1 - cor(x, y),
                        row_names_side = "right",
                        cluster_rows = F,
                        row_names_gp = gpar(fontsize = 8, col = "black"),
                        column_names_gp = gpar(fontsize = 10, col = "black"),
                        km = 6
)


library(ggpubr)
ggscatter(df, x= "Pos" , y = "Diff_G")

head(df)
table(df$Diff_G > 0)

library(tidyr)
# wide to long
gather_df <- gather(df[,c("Pos", "IDCP_G", "PCA_G", "Diff_G")], Type, Gscore, IDCP_G:Diff_G)
gather_df[1:2,]

ggscatter(gather_df, x = "Pos" , y = "Gscore", color = "Type",# alpha=0.7, 
          size = 1,
          #facet.by = "Type",
          ggtheme = theme_pubr(),
          palette = "npg")+theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))+
  xlab("Segment overlapped regions between IDCP and PCA")+
  ylab("G score")
ggsave(file.path("12.sequenza.44.samples.gistic.chrx/", "Gscore_overlapped_segment_regions_IDCP_PCA.pdf"), width = 10, height = 4)


# IDCP>PCA amplication 
amp_IDCP_big_PCA <- ol.gr.df.FDR %>% filter(IDCP_Variant_Classification == "Amp" & PCA_Variant_Classification == "Amp" & diff_gscore > 0)
dim(amp_IDCP_big_PCA)
# [1] 8  16
amp_IDCP_big_PCA





# IDCP>PCA deletion


library(karyoploteR)
kp <- plotKaryotype(genome="hg19")
kpPlotRegions(kp, ol.gr)
kpPoints(kp, chr=ol.gr.df$seqnames, x=ol.gr.df$start, y=ol.gr.df$diff_gscore,
         ymin=0, ymax=1, r0=0.05, r1=0.4, col="black", pch=".", cex=2)



# 比较基因CNV frequency
IDCP_gene <- as.data.frame(IDCP@gene.summary)
str(IDCP_gene)
head(IDCP_gene)

IDCP.gr






library(ggpubr)
df_IDCP <- as.data.frame(IDCP.gr)
head(df_IDCP)
dim(df_IDCP)
## [1] 7525   12
#df_IDCP$pos <- paste(df_IDCP$seqnames, df_IDCP$start, df_IDCP$end, sep = "_") 
df_IDCP$Cancer <- "IDCP"
table(df_IDCP$IDCP_Variant_Classification)
ggboxplot(df_IDCP, y = "IDCP_G_Score", color = "IDCP_Variant_Classification")+
  scale_color_manual(values = c("Amp" = "red", "Del" = "blue"))+
  xlab("")+
  ylab("G score (IDCP)")+
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        )

df_PCA <- as.data.frame(PCA.gr)
head(df_PCA)
dim(df_PCA)
# [1] 7111   12
df_PCA$Cancer <- "PCA"
#df_PCA$pos <- paste(df_PCA$seqnames, df_PCA$start, df_PCA$end, sep = "_")
table(PCA.gr$PCA_Variant_Classification)
ggboxplot(df_PCA, y = "PCA_G_Score", color = "PCA_Variant_Classification")+
  scale_color_manual(values = c("Amp" = "red", "Del" = "blue"))+
  xlab("")+
  ylab("G score (PCA)")+
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )



head(df_IDCP)
dim(df_IDCP)

head(df_PCA)
dim(df_PCA)

write.csv(df_IDCP, file.path("12.sequenza.44.samples.gistic.chrx/", "IDCP_gistic.csv"), quote = F, row.names = F)
write.csv(df_PCA, file.path("12.sequenza.44.samples.gistic.chrx/", "PCA_gistic.csv"), quote = F, row.names = F)



library(dplyr)
tmp <- df_PCA %>% select("PCA_Variant_Classification", "PCA_G_Score", "Cancer")
colnames(tmp) <- c("Variant_Classification", "G_Score", "Cancer")
dim(tmp)
## [1] 7111     3

tmp2 <- df_IDCP %>% select("IDCP_Variant_Classification", "IDCP_G_Score", "Cancer")
colnames(tmp2) <- c("Variant_Classification", "G_Score", "Cancer")
dim(tmp2)
## [1] 7525    3

df <- rbind(tmp, tmp2)
dim(df)
## [1] 14636     3
df$Cancer <- as.factor(df$Cancer)
df$Variant_Classification <- as.factor(df$Variant_Classification)

head(df)

# V1
ggboxplot(df, x = "Cancer", y = "G_Score", color = "Variant_Classification",
        #  add = "jitter",
         # add.params = list(size = 0.1)
          )+
  scale_color_manual(values = c("Amp" = "red", "Del" = "blue"))+
  xlab("")+
  ylab("G score")+
  theme(legend.title = element_blank() )


# V2
compare_means(G_Score ~ Cancer, data = df, group.by = "Variant_Classification")

ggboxplot(df, x = "Variant_Classification", y = "G_Score", color = "Cancer",
          #  add = "jitter",
          # add.params = list(size = 0.1)
          )+  
  stat_compare_means(aes(group = Cancer), label = "p.format", method = "t.test")+ #p.signif
  #stat_compare_means( comparisons = list(c("IDCP", "PCA")), method='wilcox.test')+
  scale_color_manual(values = c("IDCP" = "red", "PCA" = "steelblue"))+
  xlab("")+
  ylab("G score")+
  theme(legend.title = element_blank())
ggsave(file.path("12.sequenza.44.samples.gistic.chrx/", "boxplot_G_score_IDCP_PCA.pdf"), width = 4,  height = 4)


ggboxplot(df, x = "Variant_Classification", y = "G_Score", color = "Cancer",
          #  add = "jitter",
          # add.params = list(size = 0.1)
)+facet_wrap(~ Variant_Classification, scales = "free")+
  stat_compare_means(aes(group = Cancer), label = "p.format", method = "t.test")+ #p.signif
  #stat_compare_means( comparisons = list(c("IDCP", "PCA")), method='wilcox.test')+
  scale_color_manual(values = c("IDCP" = "red", "PCA" = "steelblue"))+
  xlab("")+
  ylab("G score")+
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(file.path("12.sequenza.44.samples.gistic.chrx/", "boxplot_G_score_IDCP_PCA_facet.pdf"), width = 4,  height = 4)



ggviolin(df, x = "Variant_Classification", y = "G_Score", color = "Cancer",
          #  add = "jitter",
          # add.params = list(size = 0.1)
          )+
  scale_color_manual(values = c("IDCP" = "red", "PCA" = "steelblue"))+
  xlab("")+
  ylab("G score")+
  theme(legend.title = element_blank() )




#################
kpPlotRegions(kp, losses, col="#CCFFAA")


df <- toGRanges(ol.gr.df$seqnames, ol.gr.df$start, ol.gr.df$end, ol.gr.df$diff_gscore)
df
kpPlotRegions(kp, ol.gr)


library(ggpubr)
head(ol.gr.df)
ggdensity(ol.gr.df$diff_gscore)
gghistogram(ol.gr.df, x = "diff_gscore", #add = "mean", #添加平均值线
            fill = "IDCP_fdr",
            #add_density = T,
            rug = TRUE #添加地毯线
)
            
kpSegments(kp, 
           chr=ol.gr.df$seqnames, 
           x0=ol.gr.df$start, x1=,ol.gr.df$end, 
           y0=ol.gr.df$diff_gscore-0.2, y1=ol.gr.df$diff_gscore,
           col="#AAAADD")





