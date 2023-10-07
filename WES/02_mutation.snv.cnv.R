

setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ/")
rm(list=ls())
library(maftools)
dat=read.delim2("data/46_samples.maf",sep="\t",stringsAsFactors = F,header = T)
dat[1:10,]

table(dat$Variant_Type)
dat[which(dat$Variant_Type == "MNP"), ]

tail(dat)


samples=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F)
samples$hxid_short=gsub("HXIDCP-","",samples$hxid)
samples$group[samples$group=="adenocarcinoma"]="PCA"
samples$new_id=paste(samples$hxid_short,samples$group,sep="_")
samples

pos=match(dat$Tumor_Sample_Barcode,samples$Tumor)
dat=cbind(dat,samples[pos,])
dat[1:5,]

identical(dat$sample_id, dat$Tumor)
#F
dat[which(dat$sample_id != dat$Tumor), c("sample_id", "Tumor")]

dat$Tumor_Sample_Barcode=dat$new_id

sample_order=unique(sort(dat$Tumor_Sample_Barcode))


write.table(dat,"data/46_samples.mod.maf",quote = F,sep="\t",row.names = F)



#########preprocess using maftools

maf=read.maf(maf = dat)

topgene=maf@gene.summary$Hugo_Symbol[1:20]
plotgene=unique(c(topgene,c("AR","FOXA1","NCOR2","NCOR1")))


pdf(file.path("02_data_profile","46.oncoplot.pdf"))
oncoplot(maf = maf, top = 10,showTumorSampleBarcodes = T,sampleOrder=sample_order,genes = plotgene,barcode_mar = 4,gene_mar = 6,annotationFontSize = 0.8)
dev.off()


laml.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)

pdf(file.path("02_data_profile","46.titv.pdf"))
plotTiTv(laml.titv,showBarcodes=T)
dev.off()


mutation_type=as.data.frame(laml.titv$fraction.contribution)
pos=match(mutation_type$Tumor_Sample_Barcode,samples$new_id)

mutation_type=cbind(mutation_type,samples[pos,])

library(ggpubr)

### FOR SNV
my_comparisons <- list( c("IDCP", "PCA"))
ggboxplot(mutation_type, x = "group", y = "C>T",
          add = "jitter", shape = "group",color = "group",notch = FALSE)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("C>T")


ggsave(file.path("02_data_profile","02.titv.C_to_T.pdf"),width = 6,height = 6)


my_comparisons <- list( c("IDCP", "PCA"))
ggboxplot(mutation_type, x = "group", y = "T>C",
          add = "jitter", shape = "group",color = "group",notch = FALSE)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("T>C")


ggsave(file.path("02_data_profile","02.titv.T_to_C.pdf"),width = 6,height = 6)


my_comparisons <- list( c("IDCP", "PCA"))
ggboxplot(mutation_type, x = "group", y = "C>A",
          add = "jitter", shape = "group",color = "group",notch = FALSE)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("C>A")


ggsave(file.path("02_data_profile","02.titv.C_to_A.pdf"),width = 6,height = 6)


my_comparisons <- list( c("IDCP", "PCA"))
ggboxplot(mutation_type, x = "group", y = "C>G",
          add = "jitter", shape = "group",color = "group",notch = FALSE)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("C>G")


ggsave(file.path("02_data_profile","02.titv.C_to_G.pdf"),width = 6,height = 6)


my_comparisons <- list( c("IDCP", "PCA"))
ggboxplot(mutation_type, x = "group", y = "T>A",
          add = "jitter", shape = "group",color = "group",notch = FALSE)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("T>A")


ggsave(file.path("02_data_profile","02.titv.T_to_A.pdf"),width = 6,height = 6)


my_comparisons <- list( c("IDCP", "PCA"))
ggboxplot(mutation_type, x = "group", y = "T>G",
          add = "jitter", shape = "group",color = "group",notch = FALSE)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)+ylab("T>G")


ggsave(file.path("02_data_profile","02.titv.T_to_G.pdf"),width = 6,height = 6)


pdf(file.path("02_data_profile","02.rainfallPlot.pdf"))
rainfallPlot(maf = maf, detectChangePoints = TRUE, pointSize = 1)
dev.off()




####
idcp=dat[grep("IDCP",dat$Tumor_Sample_Barcode),]
pca=dat[grep("PCA",dat$Tumor_Sample_Barcode),]

idcp_maf=read.maf(maf=idcp)
pca_maf=read.maf(maf=pca)



pt.vs.rt <- mafCompare(m1 = idcp_maf, m2 = pca_maf, m1Name = 'IDCP', m2Name = 'PCA', minMut = 2)

pdf(file.path("02_data_profile","02.forestPlot.pdf"))
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.4, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
dev.off()


pdf(file.path("02_data_profile","02.coBarplot.pdf"))
coBarplot(m1 = idcp_maf, m2 = pca_maf, m1Name = "IDCP", m2Name = "PCA")
dev.off()


pdf(file.path("02_data_profile","02.OncogenicPathways.all.pdf"))
OncogenicPathways(maf = maf)
dev.off()


pdf(file.path("02_data_profile","02.OncogenicPathways.idcp.pdf"))
OncogenicPathways(maf = idcp_maf)
dev.off()


pdf(file.path("02_data_profile","02.OncogenicPathways.pca.pdf"))
OncogenicPathways(maf = pca_maf)
dev.off()


pdf(file.path("02_data_profile","02.PlotOncogenicPathways.all.pdf"))
PlotOncogenicPathways(maf = maf, pathways = "NOTCH",showTumorSampleBarcodes = T)
dev.off()





################## CNV

cnvs=read.csv("data/0413.cnv.csv",stringsAsFactors = F)

samples_unique=unique(samples$Tumor)
genes=unique(cnvs$gene)

cnv_result=c()

for(gene in genes){
  res=c()
  for(sample in samples_unique){
    sub=cnvs[cnvs$sample_id==sample&cnvs$gene==gene,]
    if(nrow(sub)==0){
      type="NA"
    }else{
      type=sub$cnv_type
    }
    res=c(res,type)
  }
  print(length(res))
   
   cnv_result=rbind(cnv_result,res)
}

cnv_result=as.data.frame(cnv_result)
colnames(cnv_result)=samples_unique
rownames(cnv_result)=genes


pos=match(colnames(cnv_result),samples$Tumor)
colnames(cnv_result)=samples$new_id[pos]


cnv_number=apply(cnv_result,1,function(x)length(which(x!="NA")))

sorted=sort.int(cnv_number,decreasing = T,index.return = T)

cnv_df=cnv_result[sorted$ix,]


library(ComplexHeatmap)
cnv_df=as.matrix(cnv_df)

cnv_df=cnv_df[,sort(colnames(cnv_df))]

pdf(file.path("02_data_profile","02.cnv.top30.genes.pdf"),width = 10,height = 9)
Heatmap(cnv_df[1:30,], col = c("red","blue","gray"),heatmap_legend_param = list(title="CNV"))
dev.off()


#### hete:
library(mclust)

meth=c()

for(sample in samples$new_id){
  het=inferHeterogeneity(maf = maf, tsb = sample, vafCol = 'freq')
  het_mean=as.data.frame(het$clusterMeans)
  het_mean=het_mean[het_mean$meanVaf>=0.05,]
  
  tmp=c(sample,nrow(het_mean),paste(round(het_mean$meanVaf,3),collapse =";"))
  meth=rbind(meth,tmp)
}

meth=as.data.frame(meth)
colnames(meth)=c("ID","Cluster_Number","Mean_vafs")


pos=match(meth$ID,samples$new_id)
meth=cbind(meth,samples[pos,])

identical(meth$ID,meth$new_id)


pdf(file.path("02_data_profile","02.maftools.meth.example.pdf"),width = 8,height = 6)
id_example=meth$ID[which.max(meth$Cluster_Number)]
meth_example=inferHeterogeneity(maf = maf, tsb = id_example, vafCol = 'freq')
plotClusters(clusters = meth_example)
dev.off()

pdf(file.path("02_data_profile","02.maftools.meth.22_IDCP.pdf"),width = 8,height = 6)
meth_example=inferHeterogeneity(maf = maf, tsb = "22_IDCP", vafCol = 'freq')
plotClusters(clusters = meth_example)
dev.off()



meth$Cluster_Number=as.numeric(meth$Cluster_Number)

my_comparisons <- list( c("IDCP", "PCA"))
ggboxplot(meth, x = "group", y = "Cluster_Number",shape = "group",color = "group",notch = FALSE,palette = "jco",add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+ylab("Cluster Number")


ggsave(file.path("02_data_profile","02.maftools.meth.pdf"),width = 6,height = 6)



library('NMF')

mat=trinucleotideMatrix(maf)
laml.sign = estimateSignatures(mat = mat, nTry = 6)
plotCophenetic(res=laml.sign)
laml.sig = extractSignatures(mat = laml.tnm, n = 3)
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")








