rm(list=ls())

# setwd("../")
library(SomaticSignatures)
library(SomaticCancerAlterations)
#library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)


### read snvs
snvs=read.csv("data/0413.snv.csv",stringsAsFactors = F)

groups=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F)
groups$group[groups$group=="adenocarcinoma"]="PCA"

groups$hxid_group=paste(groups$hxid,groups$group,sep="_")

pos=sapply(groups$Tumor,function(x)which(snvs$sample_id==x))

dim(snvs)
snvs=snvs[unlist(pos),]
dim(snvs)

snvs=snvs[nchar(snvs$reference_allele)==1 & nchar(snvs$variant_allele)==1,]
snvs=snvs[snvs$reference_allele!="-",]
snvs=snvs[snvs$variant_allele!="-",]


dim(snvs)


pos=match(snvs$sample_id,groups$Tumor)
snvs=cbind(snvs,groups[pos,])

identical(snvs$sample_id,snvs$Tumor)


result_dir="08_mutation_signatures"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}


sca_vr = VRanges(
  seqnames =  snvs$chr ,
  ranges = IRanges(start = snvs$chr_start,end = snvs$chr_end),
  ref = snvs$reference_allele,
  alt = snvs$variant_allele,
  sampleNames = snvs$hxid_group,
  study="IDCP"
  )

sca_vr


sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19)

plotMutationSpectrum(sca_motifs, "study")
ggsave(file.path("08_mutation_signatures","plotMutationSpectrum.pdf"),width = 6,height = 3)


##
sca_mm=motifMatrix(sca_motifs, group = "sampleNames", normalize = TRUE)


### selecting number of signatures

n_sigs = 5:30
gof_nmf = assessNumberSignatures(sca_mm , n_sigs, nReplicates = 5) 

plotNumberSignatures(gof_nmf)
ggsave(file.path("08_mutation_signatures","plotNumberSignatures.pdf"),width = 6,height = 6)

save(gof_nmf,file = file.path("08_mutation_signatures","gof_nmf.Rdata"))


### select top20
sigs_nmf = identifySignatures(sca_mm,20, nmfDecomposition)

plotSignatureMap(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")
ggsave(file.path("08_mutation_signatures","plotSignatureMap.pdf"),width = 6,height = 6)


plotSignatures(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Barchart")
ggsave(file.path("08_mutation_signatures","plotSignatures.pdf"),width = 8,height = 10)

plotObservedSpectrum(sigs_nmf)
ggsave(file.path("08_mutation_signatures","plotObservedSpectrum.pdf"),width = 8,height = 15)

plotSampleMap(sigs_nmf)
ggsave(file.path("08_mutation_signatures","plotSampleMap.pdf"),width = 10,height = 10)

plotSamples(sigs_nmf)
ggsave(file.path("08_mutation_signatures","plotSamples.pdf"),width = 10,height = 8)

clu_motif = clusterSpectrum(sca_mm, "motif")

library(ggdendro)

p = ggdendrogram(clu_motif, rotate = TRUE)
p
ggsave(file.path("08_mutation_signatures","ggdendrogram.pdf"),width = 10,height = 8)


sp=sigs_nmf@signatures

colSums(sp)
sp=apply(sp,2,function(x){
  x/sum(x)
})
denovo=sp
rownames(denovo)


cosmic=read.table('https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt',header = T,sep = '\t')[,1:33]

tmp=as.character(cosmic[,2])
substr(tmp,2,2) <- '.'
rownames(cosmic)=paste(gsub('>','',cosmic[,1]),tmp)


comp=cbind(denovo[rownames(cosmic),],
           cosmic[,4:33])
colSums(comp)
pheatmap::pheatmap(cor(comp))
pheatmap::pheatmap(cor(comp)[1:20,21:50],filename = file.path("08_mutation_signatures","signature.correlation.pdf"))



