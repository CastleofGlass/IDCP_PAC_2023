
# 分开amplification和deletion比较两个IDCP和PCA
# PCA>IDCP amplification: PCA 存在的amplification cytoband，但是在IDCP中未检出
# PCA>IDCP deletion: PCA 存在的deletion cytoband，但是在IDCP中未检出
# IDCP>PCA amplification: IDCP存在的amplification cytoband，但是在PCA中未检出
# IDCP>PCA deletion: IDCP存在的deletion cytoband，但是在PCA中未检出
# 当某个amplification peak/cytoband同时在IDCP和PCA中出现时，则不考虑，认为两者间没有差异（主要缺乏重复是无法进行统计检验，或者没有阈值可以确定）；deletion peak/cytoband同理。
# 具体结果: Table Significantly different copy number alterations between IDCP and PCA
# 与第10页PPT最右热图对应，但是第10页热图包含了amplification peak/cytoband同时在IDCP和PCA中出现时类信息

# 比较IDCP和PCA各自的gistic结果文件：all_lesions.conf_90.txt
setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ")

IDCP_lesions <- read.delim("12.sequenza.44.samples.gistic.chrx.IDCP/all_lesions.conf_90.txt", header = T)
IDCP_lesions[1:2,]
dim(IDCP_lesions)
## [1] 104  32

IDCP_lesions_amp <- IDCP_lesions[1:34, c("Descriptor", "q.values", "Residual.q.values.after.removing.segments.shared.with.higher.peaks", "Wide.Peak.Limits")]
str(IDCP_lesions_amp)
colnames(IDCP_lesions_amp) <- c("Descriptor", "q.values", "Residual.q", "Wide.Peak.Limits")

IDCP_lesions_del <- IDCP_lesions[35:52, c("Descriptor", "q.values", "Residual.q.values.after.removing.segments.shared.with.higher.peaks", "Wide.Peak.Limits")]
str(IDCP_lesions_del)
colnames(IDCP_lesions_del) <- c("Descriptor", "q.values", "Residual.q", "Wide.Peak.Limits")



PCA_lesions <- read.delim("12.sequenza.44.samples.gistic.chrx.PCA/all_lesions.conf_90.txt", header = T)
PCA_lesions[1:5,]
dim(PCA_lesions)
## [1] 116  32

PCA_lesions_amp <- PCA_lesions[1:40, c("Descriptor", "q.values", "Residual.q.values.after.removing.segments.shared.with.higher.peaks", "Wide.Peak.Limits")]
str(PCA_lesions_amp)
colnames(PCA_lesions_amp) <- c("Descriptor", "q.values", "Residual.q", "Wide.Peak.Limits")

PCA_lesions_del <- PCA_lesions[41:58, c("Descriptor", "q.values", "Residual.q.values.after.removing.segments.shared.with.higher.peaks", "Wide.Peak.Limits")]
str(PCA_lesions_del)
colnames(PCA_lesions_del) <- c("Descriptor", "q.values", "Residual.q", "Wide.Peak.Limits")

head(PCA_lesions_amp)

library(dplyr)

# 1. PCA > IDCP amplification
PCA_bg_IDCP_lesions_amp <- PCA_lesions_amp %>% filter( !(Descriptor %in% IDCP_lesions_amp$Descriptor))
PCA_bg_IDCP_lesions_amp
dim(PCA_bg_IDCP_lesions_amp)
## [1] 17  4


# 2. PCA > IDCP deletion
PCA_bg_IDCP_lesions_del <- PCA_lesions_del %>% filter( !(Descriptor %in% IDCP_lesions_del$Descriptor))
PCA_bg_IDCP_lesions_del
dim(PCA_bg_IDCP_lesions_del)
## [1] 15  4

# 3. IDCP>PCA amplification （amplification: IDCP存在，但是在PCA中不在的cytoband）
IDCP_bg_PCA_lesions_amp <- IDCP_lesions_amp %>% filter( !(Descriptor %in% PCA_lesions_amp$Descriptor))
IDCP_bg_PCA_lesions_amp
dim(IDCP_bg_PCA_lesions_amp)
## [1] 9  4

# 4. IDCP>PCA deletion (deletion: IDCP存在，但是在PCA中不存在的cytoband)
IDCP_bg_PCA_lesions_del <- IDCP_lesions_del %>% filter( !(Descriptor %in% PCA_lesions_del$Descriptor))
IDCP_bg_PCA_lesions_del
dim(IDCP_bg_PCA_lesions_del)
## [1] 15  4






na_rm_PCA_bg_IDCP_lesions_amp <- na.omit(PCA_bg_IDCP_lesions_amp)
dim(na_rm_PCA_bg_IDCP_lesions_amp)
## [1] 27  7
na_rm_PCA_bg_IDCP_lesions_amp








