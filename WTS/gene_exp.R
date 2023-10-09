
library(tximport)
library(tidyverse)
library(DESeq2)


# samples infor
gene_sample <- read.table("sample_infor.xls", header = T, sep = "\t")

# read rsem result files
files <- file.path("gene_exp", paste0(gene_sample$samle_id, ".genes.results"))
names(files) <- gene_sample$samle_id
txi.rsem.g <- tximport(files, type = "rsem",txIn = FALSE, txOut = FALSE )

# prepare summary data
t_exp_anno <- data.frame(gene="",log2FoldChange="",pvalue="",sig="",group="")
diff_stat <- data.frame(group=group_infor$base,up=0,down=0)

# groups infor
group_infor <- read.table("group_infor.xls", header = T,sep = "\t")

# gene infor : ENST ENSG ENTREZID
gene_infor <- read.table("gene_infor.xls",header = T,sep = "\t")

# select genes which has entrezID
gene_infor <- gene_infor %>% dplyr::select(gene_id,entreID) %>% unique() %>% dplyr::filter(entreID != "no") 


txi.rsem.g$abundance <-
  txi.rsem.g$abundance[apply(txi.rsem.g$length,
                             1,
                             function(row) all(row !=0 )),]

txi.rsem.g$counts <-
  txi.rsem.g$counts[apply(txi.rsem.g$length,
                          1,
                          function(row) all(row !=0 )),]

txi.rsem.g$length <-
  txi.rsem.g$length[apply(txi.rsem.g$length,
                          1,
                          function(row) all(row !=0 )),]


gene_tpm <- as.data.frame(txi.rsem.g$abundance)
gene_count <- as.data.frame(txi.rsem.g$counts)

head(gene_tpm)
gene_tpm$gene_id <- rownames(gene_tpm)
gene_count$gene_id <- rownames(gene_count)


gene_tpm_anno <- merge(gene_tpm,gene_infor,by = "gene_id",all.x = T)


gene_count_anno <- merge(gene_count,gene_infor,by = "gene_id",all.x = T)


write.table(gene_tpm_anno,"gene_tpm_anno.xls",quote = F,sep = "\t",row.names = F)
write.table(gene_count_anno,"gene_count_anno.xls",quote = F,sep = "\t",row.names = F)

####################################
#------------group info---------------
sampleTable <- data.frame(condition = factor(paste0(gene_sample$sample_type)))
#   group total
sampleTable <- data.frame(condition = factor(paste0(gene_sample$sample_type,"_",gene_sample$pairPN)))
#  sub group AB 
sampleTable <- data.frame(condition = factor(paste0(gene_sample$sample_type,"_",gene_sample$class,"_",gene_sample$pairPN)))
#  sub group ADT 
sampleTable <- data.frame(condition = factor(paste0(gene_sample$sample_type,"_",gene_sample$class,"_ADT",gene_sample$ADT,"_",gene_sample$pairIP)))
sampleTable <- data.frame(condition = factor(paste0("ADT_",gene_sample$ADT,"_",gene_sample$pairIP)))


# sampleTable <- data.frame(condition = factor(gene_sample$sample_type))
rownames(sampleTable) <- colnames(txi.rsem.g$counts)



dds <- DESeqDataSetFromTximport(txi.rsem.g, sampleTable, ~condition)
dds <- DESeq(dds)
resultsNames(dds)

logfM <- data.frame(gene=rownames(results(dds)))


library(msigdbr)

H_m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

c2_m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)

c6_m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
c7_m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)


head(gene_infor)

library(tximport)
library(tidyverse)
library(DESeq2)

source('diff_exp_func.R')


############################## tumor,ctrl,prefix history list
tumor <- "IDCP_Y"
ctrl <- "NORM_Y"
prefix <- "IDCP_vs_NORM_pair"

tumor <- "PCA_Y"
ctrl <- "NORM_Y"
prefix <- "PCA_vs_NORM_pair"

tumor <- "IDCP_Y"
ctrl <- "PCA_Y"
prefix <- "IDCP_vs_PCA_pair"
#######
tumor <- "IDCP_A_Y"
ctrl <- "NORM_A_Y"
prefix <- "A_IDCP_vs_NORM_pair"

tumor <- "PCA_A_Y"
ctrl <- "NORM_A_Y"
prefix <- "A_PCA_vs_NORM_pair"

tumor <- "IDCP_A_Y"
ctrl <- "PCA_A_Y"
prefix <- "A_IDCP_vs_PCA_pair"
#######
tumor <- "IDCP_B_Y"
ctrl <- "NORM_B_Y"
prefix <- "B_IDCP_vs_NORM_pair"

tumor <- "PCA_B_Y"
ctrl <- "NORM_B_Y"
prefix <- "B_PCA_vs_NORM_pair"

tumor <- "IDCP_B_Y"
ctrl <- "PCA_B_Y"
prefix <- "B_IDCP_vs_PCA_pair"
########

tumor <- "IDCP_A_Y"
ctrl <- "IDCP_B_Y"
prefix <- "IDCP_A_vs_B_pair"

tumor <- "PCA_A_Y"
ctrl <- "PCA_B_Y"
prefix <- "PCA_A_vs_B_pair"

#############
tumor <- "IDCP_A_ADT1_Y"
ctrl <- "IDCP_A_ADT0_Y"
prefix <- "IDCP_A_ADT_1_vs_0_pair"

tumor <- "PCA_A_ADT1_Y"
ctrl <- "PCA_A_ADT0_Y"
prefix <- "PCA_A_ADT_1_vs_0_pair"

tumor <- "IDCP_B_ADT1_Y"
ctrl <- "IDCP_B_ADT0_Y"
prefix <- "IDCP_B_ADT_1_vs_0_pair"

tumor <- "PCA_B_ADT1_Y"
ctrl <- "PCA_B_ADT0_Y"
prefix <- "PCA_B_ADT_1_vs_0_pair"
#######
tumor <- "ADT_1_Y"
ctrl <- "ADT_0_Y"
prefix <- "ADT_1_vs_0_pair"

diff_exp_process(tumor,ctrl,prefix)


write.table(logfM,"group1.logFC.matrxi.xls",quote = F,sep = "\t",row.names = F)
write.table(t_exp_anno,"group1.diff_gene_exp.anno.xls",quote = F,sep = "\t",row.names = F)
write.table(diff_stat,"group1.diff_gene_exp.stat.xls",quote = F,sep = "\t",row.names = F)

head(t_exp_anno[-1,])
dif_gene <- spread(t_exp_anno[-1,c(1:2,5)],group,log2FoldChange) 
write.table(dif_gene,"group1.diff_gene_exp.xls",quote = F,sep = "\t",row.names = F)


