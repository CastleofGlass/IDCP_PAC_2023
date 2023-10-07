
# Date: 2021-11-19
# 比第一版本删除低于10%的mutation signature 
# 调整了热图的方向

# 2021-12-03
# 修改mutation signature的样本过滤标准
# HXIDCP-11的ADT信息修改

# 2021-12-09
# 更新mutation signature的结果：
# 显示signature 1 3 5 6 15 的构成比，其他的signature都合并成other
# 增加了signature的barplot, 包括proportion


rm(list=ls())

setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ/")
library(deconstructSigs)
library(magrittr)
library(openxlsx)
library(ComplexHeatmap)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(ggpubr)
library(circlize)

dat=openxlsx::read.xlsx("data/20210413_checked_21h.xlsx", sheet=1)
colnames(dat)[c(1,12,13)]=c("CP","transcript","variant_class")

### read snvs
snvs=dat[dat$variant_class=="SNV",]
snvs[1:2,]

groups=read.csv("data/12.gistic.sample.group.csv", stringsAsFactors = F) # 2021-12-03: 11号我们这周又核实了下，他是没有做ADT
groups

pos=sapply(groups$Tumor, function(x)which(snvs$sample_id==x))

dim(snvs)
## [1] 3582   13
snvs=snvs[unlist(pos),]
dim(snvs)
## [1] 2941   13

# 按照hxid对样本进行排序
snvs <- merge(snvs, groups[,c("Tumor", "hxid", "group", "ADT", "Mutation_group")], by.x = "sample_id", by.y = "Tumor")
snvs[1:6,]

snvs <- arrange(snvs, hxid)
snvs[1:6,]

snvs$group[which(snvs$group=="adenocarcinoma")] <- "PCA"

# 创建新的样本名
snvs$HXID_full <- paste(snvs$hxid, snvs$group, sep = "-")

result_dir="13_mutation_signatures_AB"
if(!file.exists(result_dir)){
  dir.create(result_dir)
}


# 使用 mut.to.sigs.input 函数，构建计算signature的输入文件，得到每个样本的96种三碱基类型。
# Convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = snvs, 
                                sample.id = "HXID_full", # "sample_id", 使用医生医院的样本名,医生要求
                                chr = "chr", 
                                pos = "chr_start", 
                                ref = "reference_allele", 
                                alt = "variant_allele")
#查看结果信息
dim(sigs.input)
#[1]  44 96  
head(t(sigs.input))

# 指定样本顺序，可以按照同一个人的PCA, IDCP样本进行排列
# samples_used=unique(snvs$sample_id)
# 修改
samples_used=unique(snvs$HXID_full)
samples_used

weights_df=c()
tumor_df=c()
product_df=c()
diff_df=c()
unknown_df=c()

# Determine the signatures contributing to samples
for(sample in samples_used){
  signatures_cosmic = whichSignatures(tumor.ref = sigs.input, 
                                      signatures.ref = signatures.cosmic, 
                                      sample.id = sample,
                                      contexts.needed = TRUE,
                                      tri.counts.method = 'default')
  weights_df=rbind(weights_df, signatures_cosmic$weights)
  tumor_df=rbind(tumor_df, signatures_cosmic$tumor)
  product_df=rbind(product_df, signatures_cosmic$product)
  diff_df=rbind(diff_df, signatures_cosmic$diff)
  unknown_df=rbind(unknown_df, c(sample, signatures_cosmic$unknown))
}

# tumor.ref：每个sample的96种三碱基突变序列
# signatures.ref：已知的signatures参考文件，可选signatures.nature2013和signatures.cosmic
# sample.id：对应tumor.ref文件中的样本名
# contexts.needed ：是否需要突变上下文
# tri.counts.method：三核酸序列标准化方式，默认“default” 不进行标准化 ；或者选择exome,genome,exome2genome,genome2exome 来限定区域。
# 通过signature.cutoff设定阈值，小于此值的为0


# whichSignatures会输出5个元素的list文件：

# weights -- data frame containing the weights assigned to each of the k signatures of the input signatures matrix
# tumor -- matrix of the trinucleotide contexts for the tumor sample used as input
# product -- matrix obtained when the tumor matrix is multiplied by the assigned weights
# diff -- matrix representing the difference between the tumor matrix and product matrix
# unknown -- numeric weight not assigned to any of the input signatures


w=lapply(samples_used , function(i){
  ## signatures.cosmic signatures.nature2013
  sample_1 = whichSignatures(tumor.ref = sigs.input[,], 
                             signatures.ref = signatures.cosmic, 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'default')
  print(i)
  return(sample_1$weights)
})

w=do.call(rbind,w)
w

library(pheatmap)
pheatmap(t(w), cluster_rows = T, cluster_cols = F)
pheatmap(w,    cluster_rows = F, cluster_cols = T)




# 指定signature权重
# 通过associated参数指定参与计算的signature
T42944S03.associate = whichSignatures(tumor.ref = sigs.input, 
                                      signatures.ref = signatures.cosmic, 
                                      sample.id = "HXIDCP-22-PCA", 
                                      associated = c("Signature.1","Signature.22"),
                                      contexts.needed = TRUE,
                                      tri.counts.method = 'default')
T42944S03.associate$weights




# plotSignatures 可视化

# Plot output
plotSignatures(T42944S03.associate, sub = 'example')
# 查看sample1的signature的组成情况，就是上面plot_example$weight , plot_example$tumor , plot_example$product 的结果可视化。

dim(weights_df) 
## [1] 44 30
weights_df[1:5,]

str(weights_df)
summary(weights_df)

# 删除在所有样本中都低于10%的signatire
# weights_df_filter <- weights_df %>% select_if(~ max(.) > 0.1)
# dim(weights_df_filter)
## [1] 44 23
rowSums(weights_df)

# 只显示signature 1 3 5 6 15 的构成比，其他的signature都合并成other
weights_df_select <- weights_df %>% dplyr::select(Signature.1, Signature.3, Signature.5, Signature.6, Signature.15)
weights_df_select[1:5,]
weights_df_noselect <- weights_df %>% dplyr::select( -c(Signature.1, Signature.3, Signature.5, Signature.6, Signature.15))
weights_df_noselect[1:5,]
weights_df_noselect$Others <- rowSums(weights_df_noselect)
weights_df_select <- cbind(weights_df_select, weights_df_noselect[,"Others"])
colnames(weights_df_select)[6] <- "Others"
rowSums(weights_df_select)

weights_df_select <- weights_df_select[order(rownames(weights_df_select)),]


df <- weights_df_select
df$HXID <- rownames(df)
df_long <- melt(df, id.vars = "HXID", variable_name = "Signatures")
head(df_long,2)
##             HXID  Signatures     value
## 1 HXIDCP-00-IDCP Signature.1 0.2443921
## 2  HXIDCP-00-PCA Signature.1 0.0000000

palette <- colorRampPalette(RColorBrewer::brewer.pal(9,name = 'Set1'))(6)
ggbarplot(df_long, x = "HXID",y = "value",
          fill = "Signatures", color = "Signatures", palette = palette,
          ylab = "Weight", xlab = ""
          #label = TRUE, lab.col = "white", lab.pos = "in"
          )+
  #scale_fill_manual(values = palette)+
  theme(axis.text.y = element_text(color = "black", size =10),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 7),
        axis.title = element_text(color = "black", size = 10))+
  ggtitle("Somatic Signatures for paired IDCP and PCA")
ggsave(file.path(result_dir, "Somatic_Signatures_for_paired_IDCP_and_PCA.pdf"), width = 6.5 , height = 4.5)


ggbarplot(df_long, x = "HXID",y = "value",
          fill = "Signatures", color = "Signatures", palette = palette,
          ylab = "Porportion", xlab = "",
          position = position_fill()
          #label = TRUE, lab.col = "white", lab.pos = "in"
)+
  #scale_fill_manual(values = palette)+
  theme(axis.text.y = element_text(color = "black", size =10),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 7),
        axis.title = element_text(color = "black", size = 10))+
  ggtitle("Somatic Signatures for paired IDCP and PCA")
ggsave(file.path(result_dir, "Somatic_Signatures_for_paired_IDCP_and_PCA_Proportion.pdf"), width = 6.5 , height = 4.5)







##################################################################################################
set.seed(13)
mat <- as.matrix(t(weights_df_select))
mat[,1:10]
colnames(mat)

dim(groups)
groups$group[which(groups$group == "adenocarcinoma")] <- "PCA"
groups$HXID_full <- paste(groups$HXID, groups$group, sep = "-")

pos <- match(colnames(mat), groups$HXID_full)
pos
df <- groups[pos,]
df
col_ha =  columnAnnotation(Cancer = df[,c('group')],
                       ADT = df[,c("ADT")],
                       Mutation_group = df[, c('Mutation_group')],
                       col=list(Mutation_group = c("A" = "#48D1CC", "B" = "#FFD700"),
                                ADT = c("0" = brewer.pal(7,"Set1")[3], "1" = brewer.pal(7,"Set1")[4]),
                                Cancer = c("IDCP" = brewer.pal(7,"Set1")[1] , "PCA" = brewer.pal(7,"Set1")[2]) ))
col_ha

f1 = colorRamp2(seq(min(mat), max(mat), length = 2), c("white", "red"))

pdf(file.path(result_dir, "Mutation_signature_mergeSignature_AB.pdf"), width = 7, height = 3.6)
Heatmap(mat, name = "weight",
        col = f1,
        cluster_rows = F,
        cluster_columns = F,
        bottom_annotation = col_ha,
        column_names_gp = gpar(fontsize = 8)
        )
dev.off()

pdf(file.path(result_dir, "Mutation_signature_mergeSignature_AB_ClusterCol.pdf"), width = 7, height = 4)
Heatmap(mat, name = "weight",
        col = f1,
        cluster_rows = F,
        cluster_columns = T,
        bottom_annotation = col_ha,
        #row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 8)
        )
dev.off()




unknown_df=as.data.frame(unknown_df, stringsAsFactors = F)
colnames(unknown_df)=c("sample","unknown")
unknown_df$unknown=as.numeric(unknown_df$unknown)

pos=match(rownames(weights_df),unknown_df$sample)

result=cbind(weights_df,unknown_df[pos,])
identical(rownames(result),result$sample)

result$sample=rownames(result)
result
dim(result)
## [1] 44 32


library(reshape)
result[1:5,]
df=melt(result, id.vars = "sample")
df[1:5,]
colnames(df)=c("sample","Signature","Weight")
pos=match(df$sample, groups$HXID_full)

df=cbind(df,groups[pos,])
df
df$group[df$group=="adenocarcinoma"]="PCA"
#write.csv(df, file.path(result_dir,"all.44samples.cosmic.signatures.csv"),quote = F,row.names = F)

dim(df)
## [1] 1364   12

ggboxplot(df, x='Signature', y= 'Weight', 
          color = "Mutation_group", add = "jitter")+
  scale_color_manual(values = c("A" = "#48D1CC", "B" = "#FFD700"))+
  #  stat_compare_means(method='anova')+
  stat_compare_means(method='t.test', aes(group = Mutation_group), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")
#ggsave(file.path(result_dir, "Boxplot_ttest_AB_group_each_signature.pdf"), width = 8, height = 4)


signatures = as.character(unique(df$Signature))
my_comparisons <- list( c("A", "B"))
pvalue_df=c()
for(cosmic_signature in signatures){
  sub=df[df$Signature==cosmic_signature,]
  pvalue=t.test(sub$Weight[sub$Mutation_group=="A"], sub$Weight[sub$Mutation_group=="B"],paired = FALSE)
  
  pvalue_df=rbind(pvalue_df,c(cosmic_signature,pvalue$p.value))
  
  if(is.na(pvalue$p.value)){
    print("NA")
  }else if(pvalue$p.value<0.05){
    pp=ggboxplot(sub,x="Mutation_group",y="Weight", shape = "Mutation_group",color = "Mutation_group",add = "point",title = cosmic_signature,ylab="Weight")+
      stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = FALSE)
    print(pp)
    ggsave(file.path(result_dir, paste("all.44samples.cosmic",cosmic_signature,"pdf",sep=".")), width = 6, height = 6)
  }
}

pvalue_df=as.data.frame(pvalue_df,stringsAsFactors = F)
colnames(pvalue_df)=c("Signature","pvalue")
pvalue_df$pvalue[pvalue_df$pvalue=="NaN"]="1"
pvalue_df$pvalue=as.numeric(pvalue_df$pvalue)

ggbarplot(pvalue_df,x="Signature",y="pvalue",sort.val = "asc",rotate=TRUE,palette = "jco",
          fill="skyblue",color = "skyblue")+ylab("T test p value between A and B")+xlab("Cosmic signatures")
#ggsave(file.path(result_dir,"all.44samples.cosmic.signatures.pvalue.between.idcp.and.pca.pdf"),width = 6,height = 6)




##################
dim(snvs)
## [1] 2941   18

pos = match(snvs$sample_id, groups$Tumor)
snvs_AB = cbind(snvs, groups[pos,])
table(snvs_AB$Mutation_group) 
##    A    B 
## 1956  985 
length(table(snvs$sample_id))
## 44

# 按照mutation group A B分组变异位点
sigs.input <- mut.to.sigs.input(mut.ref = snvs_AB, 
                                sample.id = "Mutation_group", 
                                chr = "chr", 
                                pos = "chr_start", 
                                ref = "reference_allele", 
                                alt = "variant_allele")

samples_used=unique(snvs_AB$Mutation_group)
samples_used

weights_df=c()
tumor_df=c()
product_df=c()
diff_df=c()
unknown_df=c()

for(sample in samples_used){
  signatures_cosmic = whichSignatures(tumor.ref = sigs.input, 
                                      signatures.ref = signatures.cosmic, 
                                      sample.id = sample,
                                      contexts.needed = TRUE,
                                      tri.counts.method = 'default')
  weights_df=rbind(weights_df,signatures_cosmic$weights)
  tumor_df=rbind(tumor_df,signatures_cosmic$tumor)
  product_df=rbind(product_df,signatures_cosmic$product)
  diff_df=rbind(diff_df,signatures_cosmic$diff)
  unknown_df=rbind(unknown_df,c(sample,signatures_cosmic$unknown))
}

weights_df[1:5,]

mat <- as.matrix(weights_df)
f1 = colorRamp2(seq(min(mat), max(mat), length = 2), c("white", "red"))

pdf(file.path(result_dir, "Mutation_signature_weight_real_AB_ClusterRow.pdf"), width = 5, height = 2.5)
Heatmap(mat, name = "weight",
        col = f1,
        cluster_rows = F)
dev.off()


unknown_df=as.data.frame(unknown_df,stringsAsFactors = F)
colnames(unknown_df)=c("sample","unknown")
unknown_df$unknown=as.numeric(unknown_df$unknown)

pos=match(rownames(weights_df),unknown_df$sample)

result=cbind(weights_df,unknown_df[pos,])

identical(rownames(result),result$sample)

result_df=melt(result,id.vars = "sample")
colnames(result_df)=c("Group","Signature","Weight")

ggbarplot(result_df, x = "Signature", y = "Weight",group="Group",color ="Group",position = position_dodge(0.8),palette = "npg",rotate=TRUE)
ggsave(file.path(result_dir,"Signatures.weight.between.A.and.B.pdf"),width = 4, height = 6)




