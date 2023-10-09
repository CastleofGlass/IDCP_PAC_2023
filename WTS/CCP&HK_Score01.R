
library(org.Hs.eg.db)
library(clusterProfiler)

# 31 CCP genes & 15 HK genes
# atten: No.8 gene: KIAA0101 ( Gene Synonyms: PCLAF )

CCPgenes <- c('CENPF','PTTG1','PRC1','CEP55','DLGAP5','CDCA8','ASPM','PCLAF','FOXM1','BIRC5','CDKN3','CDK1','MCM10','BUB1B','PLK1','TK1','CENPM','RAD54L','RRM2','ASF1B','KIF11','KIF20A','RAD51','CDC20','SKA1','PBK','DTL','TOP2A','NUSAP1','CDCA3','ORC6')
HKgenes <- c('RPL38','UBA52','PSMC1','RPL4','RPL37','RPS29','SLC25A3','CLTC','TXNL1','PSMA1','RPL8','MMADHC','RPL13A','PPP2CA','MRFAP1')


HKgenes_entrezID <- bitr(HKgenes,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db)

CCPgenes_entrezID <- bitr(CCPgenes,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),OrgDb = org.Hs.eg.db)

gene_tpm_anno <- read.table("gene_tpm_anno.xls",sep = "\t",header = T,stringsAsFactors = F)


