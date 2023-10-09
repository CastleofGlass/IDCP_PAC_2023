library(dplyr)
library(stringr)
library(openxlsx)
library(ComplexHeatmap)
library(ggpubr)
library(tidyverse)

# 基因列表
geneList4Score <- read.xlsx("gene_list.xlsx",sheet=1,check.names=F)
colnames(geneList4Score) <- str_remove(colnames(geneList4Score),"\\..*")

#注释信息
annoG2E <- read.table("gene2ensembl_210902.qinc",sep = "\t",header = T,stringsAsFactors = F)

# 样本对照信息
gene_sample <- read.table("sample_infor.xls", header = T, sep = "\t")
gene_sample$label <- str_c(gene_sample$huaxi,"-",gene_sample$sample_type)

#Zscore_tpm

z_score_temp <- gene_tpm_anno %>% dplyr::select(c(1:50))  %>% unique()

rownames(z_score_temp) <- z_score_temp$gene_id
z_score_mat <- z_score_temp[,-1] %>% as.matrix()


tpm_z_score <- (z_score_mat - apply(z_score_mat, 2, mean))/apply(z_score_mat, 2, sd)

tpm_z_score_df <- as.data.frame(tpm_z_score)
tpm_z_score_df$gene_id <- str_remove(rownames(tpm_z_score),"\\..*")
tpm_z_score_df <- merge(tpm_z_score_df,unique(annoG2E),by.x = "gene_id",by.y = "ensembl",all.x = T)




## XX score=====================================================================



for(i in colnames(geneList4Score)){
    
  scoreGene <- geneList4Score[,i] %>% unique()
  if(sum(is.na(scoreGene))>0){
    scoreGene <- scoreGene[-which(is.na(scoreGene))]
  }
  
  scoreGene <- scoreGene[scoreGene %in% tpm_z_score_df$symbol]
  
  score_df <- filter(tpm_z_score_df,symbol %in% scoreGene)
  rownames(score_df) <- score_df$symbol
  
  score <- apply(score_df[,-c(1,51,52)], 2, sum) %>% as.data.frame()
  score$label <- sapply(rownames(score),function(x){gene_sample$label[gene_sample$samle_id==x]})
  
  colnames(score) <- c(str_c(i,"_Score"), "Sample")
  # write.table(score,str_c(i,"_Score"),sep = "\t",col.names = T,row.names = T,quote = F)
  
  # ---------===================================================================
  
  score_tdf <- score_df[,-c(1,51,52)] %>% t()
  rownames(score_tdf) <- sapply(rownames(score_tdf),function(x){gene_sample$label[gene_sample$samle_id==x]})
  score_tdf <- as.data.frame(score_tdf)
  
  
  score_tdf$HXID <- rownames(score_tdf)
  score_tdf$pid <- str_remove(score_tdf$HXID,"-\\w+?$")
  score_tdf$type <- str_remove(score_tdf$HXID,".*-")
  score_tdf$class <- sapply(rownames(score_tdf),function(x){gene_sample$class[gene_sample$label==x]})
  score_tdf$ADT <- sapply(rownames(score_tdf),function(x){gene_sample$ADT[gene_sample$label==x]})
  score_tdf$Score <- apply(score_tdf[,1:length(scoreGene)],1,sum)
  
  write.table(score_tdf,str_c(i,"_Score.txt"),sep = "\t",row.names = T,col.names = T,quote = F)
}
  # 画组合热图===================================================================
  
  library(circlize)
  # col_fun = colorRamp2(c(-2, 0, 2), c("midnightblue", "white", "red"))
  
  
  ha = HeatmapAnnotation(class=score_tdf$class,col = list(class = c("A" = "THISTLE","B"="BISQUE")),annotation_name_side = "left")
  
  ar_count <- score_tdf[,1:length(scoreGene)] %>% t()
  test_mat <-  t(scale(t(ar_count)))

  pdf(str_c("heatmap_",i,"_Score.pdf"),width = 10,height = 8)
  Heatmap(test_mat,column_km = 3,name = "gene_exp",
          top_annotation = HeatmapAnnotation(score=anno_points(score_tdf$Score,col="red",gp = gpar(col=6)),type = score_tdf$type,class = score_tdf$class,ADT=score_tdf$ADT),
          heatmap_legend_param=list(labels_gp = gpar(fontsize = 10)),
          column_names_gp = gpar(col = c("midnightblue", "orange", "purple"), fontsize = c(10, 10,10)),
          show_row_dend = FALSE,show_parent_dend_line = FALSE,show_row_names = TRUE,row_names_gp = gpar(fontsize = 10))

  dev.off()


  #-PCA===================================================================
  library(ggplot2)
  library(grid)
  library(gridExtra)
  
  theme <- theme(panel.background = element_blank(),
                 panel.border=element_rect(fill=NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background=element_blank(),
                 axis.text.x=element_text(colour="black"),
                 axis.text.y=element_text(colour="black"),
                 axis.ticks=element_line(colour="black"),
                 plot.margin=unit(c(1,1,1,1),"line"))
  
  
  pca_temp <- score_tdf[,1:length(scoreGene)]
  rownames(pca_temp) <- score_tdf$HXID
  head(pca_temp)
  
  temp_pca <- pca_temp %>% prcomp()
  
  df_out <- as.data.frame(temp_pca$x)
  
  
  df_out$HXID <- rownames(df_out)
  df_out$type <- str_remove(df_out$HXID,".*-")
  df_out$class <- sapply(rownames(df_out),function(x){gene_sample$class[gene_sample$label==x]})
  df_out$ADT <- sapply(rownames(df_out),function(x){gene_sample$ADT[gene_sample$label==x]}) %>% as.character()
  
  # PCA percentage
  percentage <- round(temp_pca$sdev / sum(temp_pca$sdev) * 100, 2)
  percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  
  
  p0 <- ggplot(df_out,aes(x=PC1,y=PC2,color=type,shape=ADT )) +
          geom_point(size = 3) + xlab(percentage[1]) + ylab(percentage[2]) + geom_text(label = df_out$HXID,size=3,vjust = -.5,hjust=0.5)+
          theme
  ggsave(str_c("PCA_",i,"_Score0.pdf"),plot = p0,device = "pdf",width = 8,height = 8)
  
  p1 <- ggplot(df_out,aes(x=PC1,y=PC2,color=type,shape=class )) +
    geom_point(size = 3) + xlab(percentage[1]) + ylab(percentage[2]) + geom_text(label = df_out$HXID,size=3,vjust = -.5,hjust=0.5)+
    theme
  ggsave(str_c("PCA_",i,"_Score1.pdf"),plot = p1,device = "pdf",width = 8,height = 8)
  
  # feature--------------------
  df_out_r <- as.data.frame(temp_pca$rotation)
  df_out_r$feature <- as.character(row.names(df_out_r))
  
  p2 <- ggplot(df_out_r,aes(x=PC1,y=PC2,label=feature,color=feature ))+
    geom_point()+geom_text(size=3,vjust = -.5,hjust=0.5)+theme+
    theme(legend.position="none")
  
  ggsave(str_c("PCA_",i,"_Score2.pdf"),plot = p2,device = "pdf",width = 8,height = 8)
  
  
  #-箱线图--显示组间分数变化===========================================================
  library(ggpubr)
  
  
  p3 <- ggboxplot(score_tdf, x = "type", y = "Score",
            color = "type", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
            add = "jitter", shape = "type")
  my_comparisons <- list( c("PCA", "NORM"), c("IDCP", "PCA"), c("IDCP", "NORM") )
  p3 <- p3 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means() 
  
  
  p4 <- ggboxplot(score_tdf, x = "class", y = "Score",
                 color = "type", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                 add = "jitter", shape = "type") 
  
  my_comparisons <- list( c("B", "A"))
  p4 <- facet(p4 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
          stat_compare_means(), facet.by = "type")
  
  p5 <- ggboxplot(score_tdf, x = "ADT", y = "Score",
                  color = "type", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                  add = "jitter", shape = "type") 
  
  my_comparisons <- list( c("1", "0"))
  p5 <- facet(p5 + stat_compare_means(comparisons = my_comparisons,position="top")+ # Add pairwise comparisons p-value
          stat_compare_means(), facet.by = "type")
  
  ggsave(str_c("Boxplot_",i,"_Score_type.pdf"),plot = p3,device = "pdf",width = 8,height = 8)
  ggsave(str_c("Boxplot_",i,"_Score_group.pdf"),plot = p4,device = "pdf",width = 8,height = 8)
  ggsave(str_c("Boxplot_",i,"_Score_ADT.pdf"),plot = p5,device = "pdf",width = 8,height = 8)
  
  
  
  #-相关系数-显示sample和gene===========================================================
  library(corrplot)
  
  # pdf(str_c("corr_",i,"_Score1.pdf"),width = 8,height = 8)
  
  corr1 <- score_tdf[,1:length(scoreGene)] %>% t() %>% cor()
  cplot1 <- corrplot(corr1, method = "color",order="hclust", type = "upper",cl.cex= 0.8,tl.cex=0.7,addrect=4)
  
  # 按分类排列颜色
  corder <- factor(str_remove(colnames(cplot1),".*-"))
  levels(corder) <- c("#FF4500","#228B22","#1E90FF")
  
  # 出图
  pdf(str_c("corr_",i,"_Score1.pdf"),width = 8,height = 8)
  corrplot(corr1, method = "color",order="hclust", type = "upper",cl.cex= 0.8,tl.cex=0.7,addrect=4,tl.col = as.vector(corder))
  
  corr2 <- score_tdf[,1:length(scoreGene)] %>% cor()
  corrplot(corr2, method = "color",order="hclust",type = "lower", tl.col="black",cl.cex= 0.7,tl.cex=0.7,addrect=4,add=TRUE)
  
  
  corrplot(corr1, method = "color",order="hclust", type = "upper",cl.cex= 0.8,tl.cex=0.7,addrect=4)
  
  dev.off()
  
  #-折线图--显示组间分数变化===========================================================
  
  library(GGally)
  
  # pdf(str_c("parcoord_",i,"_Score.pdf"),width = 10,height = 8)
  p6 <- score_tdf %>% select(c("Score","type","pid")) %>% 
          spread(key = "type",value = "Score") %>% filter(!(is.na(IDCP)|is.na(NORM) |is.na(PCA))) %>%
            ggparcoord(columns = c(2:4),scale = "globalminmax",alphaLines = 1,groupColumn = 1,showPoints = T)  +
            scale_y_continuous()  +
            labs(x = "type",y="Score")
  # dev.off()
  ggsave(str_c("parcoord_",i,"_Score.pdf"),plot = p6,device = "pdf",width = 8,height = 8)
  
  


}

