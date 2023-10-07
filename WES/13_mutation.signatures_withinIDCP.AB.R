
rm(list=ls())

setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ/")
library(deconstructSigs)
library(magrittr)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

dat=openxlsx::read.xlsx("data/20210413_checked_21h.xlsx", sheet=1)
colnames(dat)[c(1,12,13)]=c("CP","transcript","variant_class")
dat[1:5,]

### read snvs
snvs=dat[dat$variant_class=="SNV",]
snvs[1:5,]

groups=read.csv("data/12.gistic.sample.group.csv",stringsAsFactors = F)
groups=groups[groups$group=="IDCP",]
groups[1:5,]
library(plyr)
groups <- arrange(groups, hxid)


pos=sapply(groups$Tumor, function(x)which(snvs$sample_id==x))
pos
str(pos)

dim(snvs)
## [1] 3582   13
snvs=snvs[unlist(pos),]
dim(snvs)
## [1] 1495   13


result_dir="13_mutation_signatures_withinIDCP_AB"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}

# 使用 mut.to.sigs.input 函数，构建计算signature的输入文件，得到每个样本的96种三碱基类型。
# Convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = snvs, 
                                sample.id = "sample_id", 
                                chr = "chr", 
                                pos = "chr_start", 
                                ref = "reference_allele", 
                                alt = "variant_allele")
# warning: Some samples have fewer than 50 mutations:

#查看结果信息
dim(sigs.input)
#[1]  22 96  
head(t(sigs.input))

# 指定样本顺序，可以按照同一个人的PCA, IDCP样本进行排列
samples_used=unique(snvs$sample_id)
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

# tumor.ref：每个sample的96种三碱基突变序列
# signatures.ref：已知的signatures参考文件，可选signatures.nature2013和signatures.cosmic
# sample.id：对应tumor.ref文件中的样本名
# contexts.needed ：是否需要突变上下文
# tri.counts.method：三核酸序列标准化方式，默认“default” 不进行标准化 ；或者选择exome,genome,exome2genome,genome2exome 来限定区域。

dim(weights_df) 
## [1] 22 30



######################################
weights_df[1:2,]
mat <- as.matrix(weights_df)
# 给heatmap的行做annotation, 注意对应每个变量，指定颜色
#row_ha = rowAnnotation(df=groups[,c('Mutation_group')],
#                       col=list(Mutation_group = circlize::colorRamp2(c(-1,0,1), c("blue", "white", "red"))) )

groups
row_ha = rowAnnotation(Mutation_group = groups[,c('Mutation_group')],
                       ADT = groups[,c("ADT")],
                       col=list(Mutation_group = c("A" = "#48D1CC", "B" = "#FFD700"),
                                ADT = c("0" = brewer.pal(7,"Set1")[3], "1" = brewer.pal(7,"Set1")[4])
                                ))

f1 = colorRamp2(seq(min(mat), max(mat), length = 2), c("white", "red"))

pdf(file.path(result_dir, "Mutation_signature_weight_AB.pdf"), width = 8, height = 6)
Heatmap(mat, name = "weight",
        col = f1,
        cluster_rows = F,
        right_annotation = row_ha)
dev.off()

pdf(file.path(result_dir, "Mutation_signature_weight_AB_ClusterRow.pdf"), width = 8, height = 6)
Heatmap(mat, name = "weight",
        col = f1,
        cluster_rows = T,
        right_annotation = row_ha)
dev.off()





unknown_df=as.data.frame(unknown_df, stringsAsFactors = F)
colnames(unknown_df)=c("sample","unknown")
unknown_df$unknown=as.numeric(unknown_df$unknown)

pos=match(rownames(weights_df),unknown_df$sample)
result=cbind(weights_df,unknown_df[pos,])

identical(rownames(result),result$sample)

result$sample=rownames(result)
result[1:5, ]
dim(result)
## [1] 22 32


# https://zhuanlan.zhihu.com/p/98446111

# plotSignatures 可视化
# Plot example
plot <- whichSignatures(tumor.ref = sigs.input,
                        contexts.needed = TRUE,
                        signatures.ref = signatures.cosmic,
                        sample.id = "42944S03")
samples_used

# Plot output
plotSignatures(plot, sub = '42944S03')



library(reshape)
df=melt(result, id.vars = "sample")
colnames(df)=c("sample","Signature","Weight")
df[1:5,]

pos=match(df$sample, groups$Tumor)
df=cbind(df, groups[pos,])
df$group[df$group=="adenocarcinoma"]="PCA"
df[1:5,]
dim(df)
# [1] 682  11

table(df$Signature, df$Mutation_group)
library(ggpubr)

df

ggboxplot(df, x='Signature', y= 'Weight', 
          color = "Mutation_group", add = "jitter")+
  scale_color_manual(values = c("A" = "#48D1CC", "B" = "#FFD700"))+
#  stat_compare_means(method='anova')+
  stat_compare_means(method='t.test', aes(group = Mutation_group), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")
ggsave(file.path(result_dir, "Boxplot_ttest_AB_group_each_signature.pdf"), width = 8, height = 4)


write.csv(df,file.path(result_dir,"all.22.IDCP.samples.cosmic.signatures.csv"), quote = F,row.names = F)

df[1:2,]

apply()



library(ggpubr)
signatures=as.character(unique(df$Signature))
signatures
my_comparisons <- list( c("A", "B"))
pvalue_df=c()

for(cosmic_signature in signatures){
  sub=df[df$Signature==cosmic_signature,]
  pvalue=t.test(sub$Weight[sub$Mutation_group=="A"], sub$Weight[sub$Mutation_group=="B"],paired = FALSE)
  
  pvalue_df=rbind(pvalue_df, c(cosmic_signature,pvalue$p.value))
  
  if(is.na(pvalue$p.value)){
    print("NA")
  }else if(pvalue$p.value < 0.05){
            pp=ggboxplot(sub,x="Mutation_group",y="Weight", shape = "Mutation_group",color = "Mutation_group",add = "point",title = cosmic_signature,ylab="Weight")+
            stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = FALSE)
            print(pp)
            ggsave(file.path(result_dir,paste("all.44samples.cosmic",cosmic_signature,"pdf",sep=".")),width = 6,height = 6)
  }
}

pvalue_df=as.data.frame(pvalue_df,stringsAsFactors = F)
colnames(pvalue_df)=c("Signature","pvalue")

pvalue_df$pvalue[pvalue_df$pvalue=="NaN"]="1"
pvalue_df$pvalue=as.numeric(pvalue_df$pvalue)

ggbarplot(pvalue_df,x="Signature",y="pvalue",sort.val = "asc",rotate=TRUE,palette = "jco",
          fill="skyblue",color = "skyblue")+ylab("T test p value between A and B")+xlab("Cosmic signatures")
ggsave(file.path(result_dir,"all.22.IDCP.samples.cosmic.signatures.pvalue.between.A.and.B.pdf"),width = 6,height = 6)





##################
dim(snvs)
## [1] 1495   13

pos=match(snvs$sample_id, groups$Tumor)
snvs_AB=cbind(snvs,groups[pos,])
table(snvs_AB$Mutation_group) 
length(table(snvs$sample_id))

# 按照mutation group A B分组变异位点
sigs.input <- mut.to.sigs.input(mut.ref = snvs_AB, 
                                sample.id = "Mutation_group", 
                                chr = "chr", 
                                pos = "chr_start", 
                                ref = "reference_allele", 
                                alt = "variant_allele")
samples_used=unique(snvs_AB$Mutation_group)
samples_used
## [1] "A" "B"

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
  weights_df=rbind(weights_df, signatures_cosmic$weights)
  tumor_df=rbind(tumor_df, signatures_cosmic$tumor)
  product_df=rbind(product_df, signatures_cosmic$product)
  diff_df=rbind(diff_df, signatures_cosmic$diff)
  unknown_df=rbind(unknown_df,c(sample,signatures_cosmic$unknown))
}

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
result_df
ggbarplot(result_df, x = "Signature", y = "Weight",group="Group",color ="Group",position = position_dodge(0.8),palette = "npg",rotate=TRUE)
ggsave(file.path(result_dir,"Signatures.weight.between.A.and.B.pdf"),width = 4, height = 6)






