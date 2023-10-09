
library(tximport)
library(tidyverse)
library(DESeq2)


diff_exp_process <- function(tumor,ctrl,prefix) {
  

  
  pca_res <- results(dds,contrast=c("condition",tumor,ctrl),alpha = 0.01)
  dir.create(prefix)
  outfile = paste0(prefix,"/",prefix,"_deseq2_gene.tsv")
  
  write.table(as.data.frame(pca_res),file = outfile,sep = "\t",quote = FALSE)
  diff_exp <- as.data.frame(pca_res) %>% unique() 
  
  
  
  diff_stat[diff_stat$group == prefix,]$up <- diff_exp %>% 
    dplyr::filter(pvalue < 0.01) %>% dplyr::filter(log2FoldChange > 2) %>% unique() %>%
    nrow()
  
  diff_stat[diff_stat$group == prefix,]$down <- diff_exp %>% 
    filter(pvalue < 0.01) %>% filter(log2FoldChange < -2) %>% unique() %>%
    nrow()
  
  diff_exp$gene <- rownames(diff_exp)
  logfM <- diff_exp %>% dplyr::select(gene,log2FoldChange) %>% merge(logfM,by="gene")
  colnames(logfM)[2] <- prefix
  
  
  
  diff_exp$sig <- ifelse(diff_exp$pvalue < 0.01 & diff_exp$log2FoldChange > 2,"up",
                         ifelse(diff_exp$pvalue < 0.01 & diff_exp$log2FoldChange < -2,"down","no"))
  diff_exp$group <- prefix
  temp_exp_anno <- diff_exp %>% dplyr::select(gene,log2FoldChange,pvalue,sig,group) %>% unique()
  t_exp_anno <- rbind(t_exp_anno,temp_exp_anno)
  t_exp_anno <- t_exp_anno %>% dplyr::filter(sig != "no")
  
  diff_exp_anno <- merge(diff_exp,gene_infor,by.x = "gene",by.y="gene_id") %>%
    dplyr::select(c("gene","log2FoldChange","pvalue","entreID")) %>%
    unique() 
  
  diff_exp_gene_log <- diff_exp_anno$log2FoldChange
  names(diff_exp_gene_log) <- diff_exp_anno$entreID
  diff_exp_gene_log <- sort(diff_exp_gene_log, decreasing = TRUE)
  
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  
  ggo1 <- gseGO(geneList = diff_exp_gene_log,
                OrgDb = org.Hs.eg.db,
                ont = "MF",
                pvalueCutoff = 0.05,
                verbose = TRUE)
  
  outfile = paste0(prefix,"/",prefix,"_go_mf_gsea.upset.pdf")
  pdf(outfile,width = 16,height = 10)
  upsetplot(ggo1)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_go_mf_gsea.ridge.pdf")
  pdf(outfile,width = 16,height = 10)
  ridgeplot(ggo1)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_go_mf_gsea.plot.pdf")
  pdf(outfile,width = 20,height = 10)
  gseaplot2(ggo1, geneSetID = 1:5,pvalue_table = T)
  dev.off()
  
  ggo1 <- setReadable(ggo1, org.Hs.eg.db, keyType = "ENTREZID")
  outfile = paste0(prefix,"/",prefix,"_go_mf_gsea.xls")
  write.table(ggo1,file = outfile,sep = "\t",quote = F)
  
  ggo2 <- gseGO(geneList = diff_exp_gene_log,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                verbose = TRUE)
  outfile = paste0(prefix,"/",prefix,"_go_bp_gsea.upset.pdf")
  pdf(outfile,width = 16,height = 10)
  upsetplot(ggo2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_go_bp_gsea.ridge.pdf")
  pdf(outfile,width = 16,height = 10)
  ridgeplot(ggo2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_go_mf_gsea.plot.pdf")
  pdf(outfile,width = 20,height = 10)
  gseaplot2(ggo2, geneSetID = 1:5,pvalue_table = T)
  dev.off()
  
  ggo2 <- setReadable(ggo2, org.Hs.eg.db, keyType = "ENTREZID")
  outfile = paste0(prefix,"/",prefix,"_go_bp_gsea.xls")
  write.table(ggo2,file = outfile,sep = "\t",quote = F)
  
  
  
  
  kk2 <- gseKEGG(geneList     =  diff_exp_gene_log,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  outfile = paste0(prefix,"/",prefix,"_kegg_gsea.upset.pdf")
  pdf(outfile,width = 16,height = 10)
  upsetplot(kk2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_kegg_gsea.ridge.pdf")
  pdf(outfile,width = 16,height = 10)
  ridgeplot(kk2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_kegg_gsea.plot.pdf")
  pdf(outfile,width = 20, height = 10)
  gseaplot2(kk2, geneSetID = 1:5,pvalue_table = TRUE)
  dev.off()
  
  kk2 <- setReadable(kk2, org.Hs.eg.db, keyType = "ENTREZID")
  outfile = paste0(prefix,"/",prefix,"_kegg_gsea.xls")
  write.table(kk2,file = outfile,sep = "\t",quote = F)
  
  
  
  
  
  H_em2 <- GSEA(diff_exp_gene_log, TERM2GENE = H_m_t2g)
  outfile = paste0(prefix,"/",prefix,"_msig_H_gsea.upset.pdf")
  pdf(outfile,width = 16,height = 10)
  upsetplot(H_em2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_msig_H_gsea.ridge.pdf")
  pdf(outfile,width = 16,height = 10)
  ridgeplot(H_em2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_msig_H_gsea.plot.pdf")
  pdf(outfile,width = 20, height = 10)
  gseaplot2(H_em2, geneSetID = 1:5,pvalue_table = TRUE)
  dev.off()
  
  H_em2 <- setReadable(H_em2, org.Hs.eg.db, keyType = "ENTREZID")
  outfile = paste0(prefix,"/",prefix,"_msig_H_gsea.xls")
  write.table(H_em2,file = outfile,sep = "\t",quote = F)
  
  C2_em2 <- GSEA(diff_exp_gene_log, TERM2GENE = c2_m_t2g)
  outfile = paste0(prefix,"/",prefix,"_msig_C2_gsea.upset.pdf")
  pdf(outfile,width = 16,height = 10)
  upsetplot(C2_em2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_msig_C2_gsea.ridge.pdf")
  pdf(outfile,width = 16,height = 10)
  ridgeplot(C2_em2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_msig_C2_gsea.plot.pdf")
  pdf(outfile,width = 20, height = 10)
  gseaplot2(C2_em2, geneSetID = 1:5,pvalue_table = TRUE)
  dev.off()
  
  C2_em2 <- setReadable(C2_em2, org.Hs.eg.db, keyType = "ENTREZID")
  outfile = paste0(prefix,"/",prefix,"_msig_C2_gsea.xls")
  write.table(C2_em2,file = outfile,sep = "\t",quote = F)
  
  
  
  C6_em2 <- GSEA(diff_exp_gene_log, TERM2GENE = c6_m_t2g)
  outfile = paste0(prefix,"/",prefix,"_msig_C6_gsea.upset.pdf")
  pdf(outfile,width = 16,height = 10)
  upsetplot(C6_em2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_msig_C6_gsea.ridge.pdf")
  pdf(outfile,width = 16,height = 10)
  ridgeplot(C6_em2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_msig_C6_gsea.plot.pdf")
  pdf(outfile,width = 20, height = 10)
  gseaplot2(C6_em2, geneSetID = 1:5,pvalue_table = TRUE)
  dev.off()
  
  C6_em2 <- setReadable(C6_em2, org.Hs.eg.db, keyType = "ENTREZID")
  outfile = paste0(prefix,"/",prefix,"_msig_C6_gsea.xls")
  write.table(C6_em2,file = outfile,sep = "\t",quote = F)
  
  
  C7_em2 <- GSEA(diff_exp_gene_log, TERM2GENE = c6_m_t2g)
  outfile = paste0(prefix,"/",prefix,"_msig_C7_gsea.upset.pdf")
  pdf(outfile,width = 16,height = 10)
  upsetplot(C7_em2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_msig_C7_gsea.ridge.pdf")
  pdf(outfile,width = 16,height = 10)
  ridgeplot(C7_em2)
  dev.off()
  outfile = paste0(prefix,"/",prefix,"_msig_C7_gsea.plot.pdf")
  pdf(outfile,width = 20, height = 10)
  gseaplot2(C7_em2, geneSetID = 1:5,pvalue_table = TRUE)
  dev.off()
  
  C7_em2 <- setReadable(C7_em2, org.Hs.eg.db, keyType = "ENTREZID")
  outfile = paste0(prefix,"/",prefix,"_msig_C7_gsea.xls")
  write.table(C7_em2,file = outfile,sep = "\t",quote = F)
  
  library("pathview")
  
  setwd(prefix)
  hsa05215 <- pathview(gene.data  = diff_exp_gene_log,
                       pathway.id = "hsa05215",
                       species    = "hsa",
                       out.suffix = prefix,
                       kegg.dir = "../../",
                       limit      = list(gene=6, cpd=1))
  hsa05200 <- pathview(gene.data  = diff_exp_gene_log,
                       pathway.id = "hsa05200",
                       species    = "hsa",
                       out.suffix = prefix,
                       kegg.dir = "../../",
                       limit      = list(gene=6, cpd=1))
  setwd("../")
  
  
  detach(package:clusterProfiler)
  
  
}
