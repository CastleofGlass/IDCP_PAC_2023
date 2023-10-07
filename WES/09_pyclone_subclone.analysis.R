
rm(list=ls())
options(stringsAsFactors = F)
setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ/")

files=list.files(path="09_pyclone",pattern = "loci.tsv",all.files = T,full.names = T,recursive = T)
files
result=c()

for(file in files){
  dat=read.table(file, stringsAsFactors = F, header = T, sep="\t")
  sample_ids = unique(dat$sample_id)
  
  cluster_vars = table(dat$cluster_id)
  clusters = cluster_vars[cluster_vars >= 3*length(sample_ids)]  ## at least 3 variants in the clusters
  
  pos=sapply(names(clusters), function(x){which(dat$cluster_id==x)})
  
  print(nrow(dat))
  dat = dat[unlist(pos), ]
  print(nrow(dat))
  
  for(id in sample_ids){
    sub=dat[dat$sample_id == id,]
    
    # 计算每个cluster的ccf均值，过滤ccf低于0.05的cluster 
    mean_ccf = tapply(sub$cellular_prevalence, sub$cluster_id, mean)
    mean_ccf = mean_ccf[mean_ccf > 0.05]
    
    # 获取保留下的cluster在原来数据中行的位置 
    sub_pos = sapply(names(mean_ccf), function(x) which(sub$cluster_id == x) ) 
    sub = sub[unlist(sub_pos), ]
    
    sub$clonalty = "Subclone"
    
    # mean ccf最大的cluster 即为主克隆
    clone_number = names(mean_ccf)[which.max(mean_ccf)]
    sub$clonalty[sub$cluster_id == clone_number] = "Clone"
    
    result=rbind(result,sub)
  }
}
result = as.data.frame(result)
result[1:5,]
result$gene = sapply(result$mutation_id, function(x) unlist(strsplit(x, ":"))[5] )

samples = read.csv("data/idcp.samples.with.HXid.update.csv", stringsAsFactors = F)
samples

pos = sapply(samples$Tumor, function(x) which(result$sample_id == x))
result = result[unlist(pos), ]
pos = match(result$sample_id, samples$Tumor)
result=cbind(result, samples[pos,])
identical(result$sample_id, result$Tumor)
result$group[result$group == "adenocarcinoma"]="PCA"
result$hxid_group = paste(result$hxid, result$group, sep="_")
result[1:5, ]

write.csv(result,file.path("09_pyclone","clone_subclone_classification.csv"),quote = F,row.names = F)

hxids = unique(result$hxid_group)
genes = unique(result$gene)

heatmap_df=c()
for(gene in genes){
  gene_res=c()
  for(id in hxids){
    tmp = result[result$hxid_group==id & result$gene==gene, ]
    if(nrow(tmp)==0){
      type=""
    }else{
      type=paste(unique(sort(tmp$clonalty)), collapse = "_")
    }
    gene_res=c(gene_res, type)
  }
  heatmap_df=rbind(heatmap_df,gene_res)
}

#heatmap_df=as.data.frame(heatmap_df)
rownames(heatmap_df)=genes
colnames(heatmap_df)=hxids

number = apply(heatmap_df, 1 ,function(x)length(which(x=="")))
ord = sort.int(number, decreasing = F, index.return = T)
ord$ix
heatmap_df = heatmap_df[ord$ix,]
write.csv(heatmap_df,file.path("09_pyclone","clone_subcone_summary.csv"),row.names = T,quote = F)

library(ComplexHeatmap)
library(RColorBrewer)
library(ggpubr)
pdf(file.path("09_pyclone","clone_subclone_top30.gene.pdf"),width = 10,height = 7)
Heatmap(heatmap_df[1:30,],
        heatmap_legend_param=list(title=""),
        col=brewer.pal(4,"Set1"))
dev.off()


num_df=c()

for(i in 1:30){
  clone_num=length(which(heatmap_df[i,] == "Clone"))
  subclone_num=length(which(heatmap_df[i,] == "Subclone"))
  clone_subclone=length(which(heatmap_df[i,] == "Clone_Subclone"))
  num_df=rbind(num_df, c(rownames(heatmap_df)[i], clone_num,subclone_num, clone_subclone))
}

num_df = as.data.frame(num_df, stringsAsFactors = F)
colnames(num_df)=c("Gene","Clone","Subclone","Clone_sublone")

num_df[,2:4] = apply(num_df[,2:4], 2, as.numeric)

library(reshape)
df = melt(num_df, id.vars = "Gene")
colnames(df) = c("Gene", "Type", "Number")

ggbarplot(df, "Gene", "Number",
                fill = "Type", color = "Type",palette="jco",order = rev(num_df$Gene),
                label = FALSE)+rotate()+ylab("Variants Number")
ggsave(file.path("09_pyclone","clone_subclone.gene.number.barplot.pdf"))


### clone subclone sigantures
result$chr=sapply(result$mutation_id, function(x)unlist(strsplit(x,":"))[1])
result$start=as.numeric(sapply(result$mutation_id,function(x)unlist(strsplit(x,":"))[2]))
result$ref=sapply(result$mutation_id,function(x)unlist(strsplit(x,":"))[3])
result$alt=sapply(result$mutation_id,function(x)unlist(strsplit(x,":"))[4])


#remove indel
result=result[nchar(result$ref)==1&nchar(result$alt)==1,]
result=result[result$ref!="-",]
result=result[result$alt!="-",]

library(deconstructSigs)


sigs.input <- mut.to.sigs.input(mut.ref = result, 
                                sample.id = "clonalty", 
                                chr = "chr", 
                                pos = "start", 
                                ref = "ref", 
                                alt = "alt")

samples_used=unique(result$clonalty)

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
  pdf(file.path("09_pyclone",paste(sample,"signautre.pie.pdf",sep=".")))
  makePie(signatures_cosmic)
  dev.off()
  
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


signautre_df=cbind(weights_df,unknown_df[pos,])

identical(rownames(signautre_df),signautre_df$sample)

signautre_df$sample=rownames(signautre_df)

df=melt(signautre_df,id.vars = "sample")

colnames(df)=c("sample","Signature","Weight")

p=ggbarplot(df, x="sample", y="Weight",color = "Signature",fill="Signature")+rotate_x_text(angle = 90)

pdf(file.path("09_pyclone","clone_subclone_signautre.percentage.barplot.pdf"))
print(p)
dev.off()












result$newid=paste(result$group,result$clonalty,sep="_")



library(deconstructSigs)


sigs.input <- mut.to.sigs.input(mut.ref = result, 
                                sample.id = "newid", 
                                chr = "chr", 
                                pos = "start", 
                                ref = "ref", 
                                alt = "alt")




samples_used=unique(result$newid)

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
  pdf(file.path("09_pyclone",paste(sample,"signautre.pie.pdf",sep=".")))
  makePie(signatures_cosmic)
  dev.off()
  
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


signautre_df=cbind(weights_df,unknown_df[pos,])

identical(rownames(signautre_df),signautre_df$sample)

signautre_df$sample=rownames(signautre_df)

df=melt(signautre_df,id.vars = "sample")

colnames(df)=c("sample","Signature","Weight")

p=ggbarplot(df, x="sample", y="Weight",color = "Signature",fill="Signature")+rotate_x_text(angle = 90)

pdf(file.path("09_pyclone","IDCP_vs_PCA.clonality.signautre.percentage.barplot.pdf"))
print(p)
dev.off()


