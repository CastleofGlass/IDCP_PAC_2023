
rm(list=ls())
options(stringsAsFactors = F)

# setwd("C:\\Users\\Yunjie.song\\Desktop\\20210408_IDCP")
setwd("/Volumes/Temp/医学肿瘤研究套路/IDCP/SYJ")

samples=read.csv("data/idcp.pairs.with.HXid.csv",stringsAsFactors = F)
samples
library(clonevol)

driver_genes=read.table("mobster_378_driver_genes_from_Martincorena_and_Tarabichi.txt",stringsAsFactors = F,header = F)
driver_genes=driver_genes$V1
driver_genes

args="HXIDCP-18"

normal=samples$normal[samples$hxid==args]

files=list.files(path=file.path("10_72425_pyclone",normal,"tables"),pattern = "loci.tsv",all.files = T,full.names = T)
files
## [1] "10_72425_pyclone/72425S03X2/tables/loci.tsv"

dat=read.table(files,stringsAsFactors = F,sep="\t",header = T)
dat[1:3,]
table(dat$sample_id)
## 72425S02X2   72425S04   72425S08   72425S09 
##         69         69         69         69 
dim(dat)

samples_used=unique(dat$sample_id)
samples_used

cluster_numbers=table(dat$cluster_id)
cluster_numbers
##   0   1   2   3   4   5   6   7 
## 140  32  76   8   4   8   4   4 

#cluster_numbers=cluster_numbers[cluster_numbers>=3*length(samples_used)]  # only used clusters with at least 3 variants
#cluster_numbers

# 2021-12-23
cluster_numbers=cluster_numbers[cluster_numbers>=4]  # only used clusters with at least 4 variants
cluster_numbers

pos=sapply(names(cluster_numbers),function(x)which(dat$cluster_id==x))

dat_used=dat[unlist(pos),]
dat_used

samples_df=read.csv("data/idcp.samples.with.HXid.72425.csv",stringsAsFactors = F)
samples_df
##       normal      Tumor    id      hxid          group
## 1 72425S03X2 72425S02X2 72425 HXIDCP-18           IDCP
## 2 72425S03X2   72425S04 72425 HXIDCP-18 adenocarcinoma
## 3 72425S03X2   72425S08 72425 HXIDCP-18         lymph1
## 4 72425S03X2   72425S09 72425 HXIDCP-18         lymph2


## change name
pos2=match(dat_used$sample_id,samples_df$Tumor)
dat_used=cbind(dat_used, samples_df[pos2,])
dat_used

dat_used$gene=sapply(dat_used$mutation_id, function(x) unlist(strsplit(x,":"))[5])

dat_used$driver_status=sapply(dat_used$gene,function(x){if(x %in% driver_genes){return(TRUE)}else{return(FALSE)}})

dat_used_clusters_mean=tapply(dat_used$cellular_prevalence,dat_used$cluster_id, mean)
dat_used_clusters_mean
dat_used_clusters_mean_sorted=sort(dat_used_clusters_mean, decreasing = T)
dat_used_clusters_mean_sorted

id_dict=seq(1:length(dat_used_clusters_mean_sorted))
id_dict
names(id_dict)=names(dat_used_clusters_mean_sorted)
id_dict
## 1 0 2 
## 1 2 3 

## 3 5 6 4 7 1 0 2 
## 1 2 3 4 5 6 7 8 

# 2021-12-23
## 3 5 6 4 7 1 0 2 
## 1 2 3 4 5 6 7 8 


dat_used$new_cluster = id_dict[as.character(dat_used$cluster_id)]
dat_used[1:3,]
dat_used$group[dat_used$group=="adenocarcinoma"]="PCA"

result=c()
variants=unique(dat_used$mutation_id)

for(variant in variants){
  sub=dat_used[dat_used$mutation_id == variant,]
  cluster=unique(sub$new_cluster)
  gene=unique(sub$gene)
  is.driver=unique(sub$driver_status)
  groups=sort(unique(sub$group))
  
  ccfs=c()
  vafs=c()
  for(abb in groups){
    ccf=100*sub$cellular_prevalence[sub$group==abb]
    ccfs=c(ccfs,ccf)
    
    vaf=100*sub$variant_allele_frequency[sub$group==abb]
    vafs=c(vafs,vaf)
  }
  
  names(ccfs)=paste(groups,"ccf",sep=".")
  names(vafs)=paste(groups,"vaf",sep=".")
  
  res=c(cluster, gene, is.driver, variant, ccfs,vafs)
  names(res)[1:4] = c("cluster", "gene", "is.driver", "variant_id")
  result=rbind(result, res) 
}

result=as.data.frame(result, stringsAsFactors = F)
result[1:4,]
dim(result)
## [1] 62 12

vaf_col_name = grep('.vaf', colnames(result), value=T)
ccf_col_name = grep('.ccf',colnames(result), value = T)

result$cluster = as.numeric(result$cluster)
result[,5:ncol(result)] = apply(result[,5:ncol(result)], 2, as.numeric)
result$is.driver = as.logical(result$is.driver)


### use ccf

pdf(paste("10_72425_pyclone/",args,".variant.cluster.ccf.20211223.pdf",sep=""), width = 7, height = 7, useDingbats = FALSE, title='')
pp <- plot.variant.clusters(result,
                            cluster.col.name = 'cluster',
                            show.cluster.size = FALSE,
                            cluster.size.text.color = 'blue',
                            vaf.col.names = ccf_col_name)
dev.off()



### use vaf

pdf(paste("10_72425_pyclone/",args,".variant.cluster.vaf.20211223.pdf",sep=""), width = 7, height = 7, useDingbats = FALSE, title='')
pp <- plot.variant.clusters(result,
                            cluster.col.name = 'cluster',
                            show.cluster.size = FALSE,
                            cluster.size.text.color = 'blue',
                            vaf.col.names = vaf_col_name)
dev.off()


#### use ccf
pdf(paste("10_72425_pyclone/",args,".cluster.flow.ccf.20211223.pdf",sep=""),width = 7,height = 5)
plot.cluster.flow(result, vaf.col.names = ccf_col_name,
                  sample.names = ccf_col_name)
dev.off()

#### use vaf
pdf(paste("10_72425_pyclone/",args,".cluster.flow.vaf.20211223.pdf",sep=""),width = 7,height = 5)
plot.cluster.flow(result, vaf.col.names = vaf_col_name,
                  sample.names = vaf_col_name)
dev.off()


plot.pairwise(result, col.names = vaf_col_name,
              out.prefix = paste("10_72425_pyclone/",args,".variants.pairwise.plot.vaf",sep="")
              )

plot.pairwise(result, col.names = ccf_col_name,
              out.prefix = paste("10_72425_pyclone/",args,".variants.pairwise.plot.ccf.20211223",sep="")
)




####### mono vaf
y_vaf = infer.clonal.models(variants = result,
                        cluster.col.name = 'cluster',
                        vaf.col.names =  vaf_col_name ,
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = 1,
                        cluster.center = 'median',
                        ignore.clusters = NULL,
                        min.cluster.vaf = 0.05,
                        # min probability that CCF(clone) is non-negative
                        sum.p = 0.05,
                        # alpha level in confidence interval estimate for CCF(clone)
                        alpha = 0.05
)


y_vaf <- transfer.events.to.consensus.trees(y_vaf,
                                        result[result$is.driver, ],
                                        cluster.col.name = 'cluster',
                                        event.col.name = 'gene')

y_vaf <- convert.consensus.tree.clone.to.branch(y_vaf, branch.scale = 'sqrt',cluster.col = "cluster")


plot.clonal.models(y_vaf,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 0.0001,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   max.num.models.to.plot=10,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = file.path("10_72425_pyclone/",paste(args,"mono_vaf_20211223",sep="_")),
                   out.format = 'pdf',
                   overwrite.output = TRUE,width = 12,
                   height = 6,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,6))









####### mono ccf
y_ccf = infer.clonal.models(variants = result,
                            cluster.col.name = 'cluster',
                            vaf.col.names =  ccf_col_name ,
                            cancer.initiation.model='monoclonal',
                            subclonal.test = 'bootstrap',
                            subclonal.test.model = 'non-parametric',
                            num.boots = 1000,
                            founding.cluster = 1,
                            cluster.center = 'median',
                            ignore.clusters = NULL,
                            min.cluster.vaf = 0,
                            # min probability that CCF(clone) is non-negative
                            sum.p = 0.05,
                            # alpha level in confidence interval estimate for CCF(clone)
                            alpha = 0.05
)


y_ccf <- transfer.events.to.consensus.trees(y_ccf,
                                            result[result$is.driver,],
                                            cluster.col.name = 'cluster',
                                            event.col.name = 'gene')

y_ccf <- convert.consensus.tree.clone.to.branch(y_ccf, branch.scale = 'sqrt',cluster.col = "cluster")


plot.clonal.models(y_ccf,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 0.0001,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   max.num.models.to.plot=10,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = file.path("10_72425_pyclone",paste(args,"mono_ccf_20211223",sep="_")),
                   out.format = 'pdf',
                   overwrite.output = TRUE,width = 12,
                   height = 6,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,6))






####### poly vaf
y_vaf = infer.clonal.models(variants = result,
                            cluster.col.name = 'cluster',
                            vaf.col.names =  vaf_col_name ,
                            cancer.initiation.model='polyclonal',
                            subclonal.test = 'bootstrap',
                            subclonal.test.model = 'non-parametric',
                            num.boots = 1000,
                            founding.cluster = 1,
                            cluster.center = 'median',
                            ignore.clusters = NULL,
                            min.cluster.vaf = 0.05,
                            # min probability that CCF(clone) is non-negative
                            sum.p = 0.05,
                            # alpha level in confidence interval estimate for CCF(clone)
                            alpha = 0.05
)

y_vaf <- transfer.events.to.consensus.trees(y_vaf,
                                            result[result$is.driver,],
                                            cluster.col.name = 'cluster',
                                            event.col.name = 'gene')

y_vaf <- convert.consensus.tree.clone.to.branch(y_vaf, branch.scale = 'sqrt',cluster.col = "cluster")

plot.clonal.models(y_vaf,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 0.0001,
                   mtcab.show.event=TRUE,
                   mtcab.tree.rotation.angle=180,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   max.num.models.to.plot=10,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = file.path("10_72425_pyclone/",paste(args,"poly_vaf_20211223",sep="_")),
                   out.format = 'pdf',
                   overwrite.output = TRUE,width = 12,
                   height = 6,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,8))









####### poly ccf
#pos1=which(result$LAM.vaf>50)
#pos2=which(result$LBM.vaf>50)
#pos3=which(result$SCP.vaf>50)

#result_poly_ccf=result[-unique(c(pos1,pos2,pos3)),]

y_ccf = infer.clonal.models(variants = result,
                            cluster.col.name = 'cluster',
                            vaf.col.names =  ccf_col_name ,
                            cancer.initiation.model='polyclonal',
                            subclonal.test = 'bootstrap',
                            subclonal.test.model = 'non-parametric',
                            num.boots = 1000,
                            founding.cluster = 1,
                            cluster.center = 'median',
                            ignore.clusters = NULL,
                            min.cluster.vaf = 0.05,
                            # min probability that CCF(clone) is non-negative
                            sum.p = 0.05,
                            # alpha level in confidence interval estimate for CCF(clone)
                            alpha = 0.05
)


y_ccf <- transfer.events.to.consensus.trees(y_ccf,
                                            result[result$is.driver,],
                                            cluster.col.name = 'cluster',
                                            event.col.name = 'gene')

y_ccf <- convert.consensus.tree.clone.to.branch(y_ccf, branch.scale = 'sqrt',cluster.col = "cluster")


plot.clonal.models(y_ccf,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 0.0001,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   max.num.models.to.plot=10,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = file.path("10_72425_pyclone/",paste(args,"poly_ccf_20211223",sep="_")),
                   out.format = 'pdf',
                   overwrite.output = TRUE,width = 12,
                   height = 6,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,6))


