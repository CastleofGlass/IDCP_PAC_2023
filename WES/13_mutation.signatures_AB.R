

rm(list=ls())

# setwd("../")
library(deconstructSigs)
library(magrittr)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

dat=openxlsx::read.xlsx("data/20210413_checked_21h.xlsx", sheet=1)
colnames(dat)[c(1,12,13)]=c("CP","transcript","variant_class")

### read snvs
snvs=dat[dat$variant_class=="SNV",]

groups=read.csv("data/12.gistic.sample.group.csv",stringsAsFactors = F)
groups

pos=sapply(groups$Tumor,function(x)which(snvs$sample_id==x))

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

library(plyr)

result_dir="13_mutation_signatures_AB"

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
#查看结果信息
dim(sigs.input)
#[1]  44 96  
head(t(sigs.input))

# 指定样本顺序，可以按照同一个人的PCA, IDCP样本进行排列
samples_used=unique(snvs$sample_id)
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
library(pheatmap)
pheatmap(t(w),cluster_rows = T, cluster_cols = F)
pheatmap(w,   cluster_rows = F, cluster_cols = T)


# 指定signature权重
# 通过associated参数指定参与计算的signature
T42944S03.associate = whichSignatures(tumor.ref = sigs.input, 
                                     signatures.ref = signatures.cosmic, 
                                     sample.id = "42944S03", 
                                     associated = c("Signature.1","Signature.22"),
                                     contexts.needed = TRUE,
                                     tri.counts.method = 'default')
T42944S03.associate$weights


# plotSignatures 可视化

# Plot output
plotSignatures(T42944S03.associate, sub = 'example')
# 查看sample1的signature的组成情况，就是上面plot_example$weight , plot_example$tumor , plot_example$product 的结果可视化。



weights_df[1:5,]

set.seed(13)
mat <- as.matrix(weights_df)

pos <- match(rownames(mat), groups$Tumor)
df <- groups[pos,]

f1 = colorRamp2(seq(min(mat), max(mat), length = 2), c("white", "red"))
row_ha = rowAnnotation(Cancer = df[,c('group')],
                       ADT = df[,c("ADT")],
                       Mutation_group = df[, c('Mutation_group')],
                       col=list(Mutation_group = c("A" = "#48D1CC", "B" = "#FFD700"),
                                ADT = c("0" = brewer.pal(7,"Set1")[3], "1" = brewer.pal(7,"Set1")[4]),
                                Cancer = c("IDCP" = brewer.pal(7,"Set1")[1] , "adenocarcinoma" = brewer.pal(7,"Set1")[2]) ))
row_ha
pdf(file.path(result_dir, "Mutation_signature_weight_AB.pdf"), width = 8, height = 8)
Heatmap(mat, name = "weight",
        col = f1,
        cluster_rows = F,
        right_annotation = row_ha)
dev.off()




unknown_df=as.data.frame(unknown_df,stringsAsFactors = F)
colnames(unknown_df)=c("sample","unknown")
unknown_df$unknown=as.numeric(unknown_df$unknown)


pos=match(rownames(weights_df),unknown_df$sample)


result=cbind(weights_df,unknown_df[pos,])

identical(rownames(result),result$sample)

result$sample=rownames(result)


library(reshape)


df=melt(result,id.vars = "sample")

colnames(df)=c("sample","Signature","Weight")

pos=match(df$sample,groups$Tumor)

df=cbind(df,groups[pos,])
df$group[df$group=="adenocarcinoma"]="PCA"


write.csv(df,file.path(result_dir,"all.44samples.cosmic.signatures.csv"),quote = F,row.names = F)







library(ggpubr)

signatures=as.character(unique(df$Signature))

my_comparisons <- list( c("A", "B"))

pvalue_df=c()

for(cosmic_signature in signatures){
  sub=df[df$Signature==cosmic_signature,]
  pvalue=t.test(sub$Weight[sub$Mutation_group=="A"],sub$Weight[sub$Mutation_group=="B"],paired = FALSE)
  
  pvalue_df=rbind(pvalue_df,c(cosmic_signature,pvalue$p.value))
  
  if(is.na(pvalue$p.value)){
    print("NA")
  }else if(pvalue$p.value<0.05){
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

ggsave(file.path(result_dir,"all.44samples.cosmic.signatures.pvalue.between.idcp.and.pca.pdf"),width = 6,height = 6)







##################
dim(snvs)

pos=match(snvs$sample_id,groups$Tumor)

snvs_AB=cbind(snvs,groups[pos,])


sigs.input <- mut.to.sigs.input(mut.ref = snvs_AB, 
                                sample.id = "Mutation_group", 
                                chr = "chr", 
                                pos = "chr_start", 
                                ref = "reference_allele", 
                                alt = "variant_allele")


samples_used=unique(snvs_AB$Mutation_group)

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


unknown_df=as.data.frame(unknown_df,stringsAsFactors = F)
colnames(unknown_df)=c("sample","unknown")
unknown_df$unknown=as.numeric(unknown_df$unknown)


pos=match(rownames(weights_df),unknown_df$sample)


result=cbind(weights_df,unknown_df[pos,])

identical(rownames(result),result$sample)




result_df=melt(result,id.vars = "sample")
colnames(result_df)=c("Group","Signature","Weight")

ggbarplot(result_df,"Signature","Weight",group="Group",color ="Group",position = position_dodge(0.8),palette = "jco",rotate=TRUE)
ggsave(file.path(result_dir,"Signatures.weight.between.A.and.B.pdf"),width = 6,height = 6)


