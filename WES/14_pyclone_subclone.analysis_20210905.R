rm(list=ls())
options(stringsAsFactors = F)

##  setwd("../")
result_dir="14_subclone_reanalysis_202109"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}


groups=read.csv("data/12.gistic.sample.group.csv",stringsAsFactors = F)
groups$group[groups$group=="adenocarcinoma"]="PCA"

groups$id_pair=paste(groups$HXID,groups$group,sep="_")



df=read.csv("09_pyclone/clone_subcone_summary.csv",stringsAsFactors = F,check.names = F,row.names = 1)


############## fisher exact test for IDCP

mut_sample <- groups[groups$group == "IDCP", "id_pair"]
wild_sample <- groups[groups$group == "PCA", "id_pair"]



#### function 
fisherresult=function(focal_genes,mut_sample,wild_sample){
  
  mut=as.matrix(focal_genes[,mut_sample])
  wild=as.matrix(focal_genes[,wild_sample])
  
  
  library(future.apply)
  plan(multiprocess)
  p_values <- future_lapply(seq(nrow(mut)), function(x){
    
    #res <- wilcox.test(x = mut[x,], y =  wild[x,],paired = isPaired)
    mut_snv=mut[x,]
    mut_clone=length(grep("Clone",mut_snv))
    mut_subclone=length(grep("Subclone",mut_snv))
    
    wild_snv=wild[x,]
    
    wild_clone=length(grep("Clone",wild_snv))
    wild_subclone=length(grep("Subclone",wild_snv))
    
    mat=cbind(c(mut_clone,mut_subclone),c(wild_clone,wild_subclone))
    colnames(mat)=c("Mut","Wild")
    rownames(mat)=c("Clone","Sublone")
    
    res=fisher.test(mat)
    return(res$p.value)
    
  })
  
  
  p <- unlist(p_values)
  padj=p.adjust(p,method = "fdr")
  
  
  df <- data.frame(Symbol = rownames(mut),
                   #Mut_clone=mat[1,1],
                   #Mut_Subclone=mat[2,1],
                   #Wild_Clone=mat[2,1],
                   #Wild_Subclone=mat[2,2],
                   pvalue = p,
                   padj=padj,
                   logp=-log10(p),stringsAsFactors = F)
  
  
  df$significance="No"
  df$significance[df$pvalue<0.05]="Yes"
  
  return(df)
}


idcp_pca_clone_subclone=fisherresult(df,mut_sample,wild_sample)

write.csv(idcp_pca_clone_subclone,file.path(result_dir,"IDCP_PCA_clone_subclone.csv"),quote = F,row.names = F)






############## fisher exact test for ADT

mut_sample <- groups[groups$ADT == "1", "id_pair"]
wild_sample <- groups[groups$ADT == "0", "id_pair"]


ADT_clone_subclone=fisherresult(df,mut_sample,wild_sample)

write.csv(ADT_clone_subclone,file.path(result_dir,"ADT_clone_subclone.csv"),quote = F,row.names = F)









############## fisher exact test for group A B

mut_sample <- groups[groups$Mutation_group == "A", "id_pair"]
wild_sample <- groups[groups$Mutation_group == "B", "id_pair"]


Group_AB_clone_subclone=fisherresult(df,mut_sample,wild_sample)

write.csv(Group_AB_clone_subclone,file.path(result_dir,"Group_AB_clone_subclone.csv"),quote = F,row.names = F)




#### with in IDCP  A_VS_B

groups_idcp=groups[groups$group=="IDCP",]

mut_sample <- groups_idcp[groups_idcp$Mutation_group == "A", "id_pair"]
wild_sample <- groups_idcp[groups_idcp$Mutation_group == "B", "id_pair"]


IDCP_A_vs_B_clone_subclone=fisherresult(df,mut_sample,wild_sample)
write.csv(IDCP_A_vs_B_clone_subclone,file.path(result_dir,"IDCP_A_vs_B_clone_subclone.csv"),quote = F,row.names = F)






#### with in PCA  A_VS_B

groups_pca=groups[groups$group=="PCA",]

mut_sample <- groups_pca[groups_pca$Mutation_group == "A", "id_pair"]
wild_sample <- groups_pca[groups_pca$Mutation_group == "B", "id_pair"]


PCA_A_vs_B_clone_subclone=fisherresult(df,mut_sample,wild_sample)
write.csv(PCA_A_vs_B_clone_subclone,file.path(result_dir,"PCA_A_vs_B_clone_subclone.csv"),quote = F,row.names = F)



#### row annotation

gene1=Group_AB_clone_subclone[sort.int(Group_AB_clone_subclone$pvalue,decreasing = F,index.return = T)$ix,]
gene2=IDCP_A_vs_B_clone_subclone[sort.int(IDCP_A_vs_B_clone_subclone$pvalue,decreasing = F,index.return = T)$ix,]
gene3=PCA_A_vs_B_clone_subclone[sort.int(PCA_A_vs_B_clone_subclone$pvalue,decreasing = F,index.return = T)$ix,]


top_10_genes=unique(c(gene1$Symbol[1:10],gene2$Symbol[1:10],gene3$Symbol[1:10]))

anno_df=data.frame(Gene=top_10_genes,stringsAsFactors = F)
anno_df$A_vs_B=gene1$significance[match(top_10_genes,gene1$Symbol)]
anno_df$IDCP_AB=gene2$significance[match(top_10_genes,gene2$Symbol)]
anno_df$PCA_AB=gene3$significance[match(top_10_genes,gene3$Symbol)]


rownames(anno_df)=anno_df$Gene
anno_df$Gene=NULL



### col annotation

col_anno=data.frame(Name=colnames(df),stringsAsFactors = F)

groups_update=groups[match(col_anno$Name,groups$id_pair),]

identical(col_anno$Name,groups_update$id_pair)

col_anno$ADT=as.character(groups_update$ADT)
col_anno$AB=groups_update$Mutation_group
rownames(col_anno)=col_anno$Name
col_anno$Name=NULL
colnames(col_anno)=paste("Group",colnames(col_anno),sep="_")

df_plot=df[top_10_genes,]


library(ComplexHeatmap)
library(RColorBrewer)
library(ggpubr)



pdf(file.path(result_dir,"top_gene_betwwen_A_vs_B.pdf"),width = 12,height = 7)
Heatmap(as.matrix(df_plot),heatmap_legend_param=list(title=""),col=brewer.pal(4,"Set1"),
        left_annotation = rowAnnotation(df=anno_df),
        top_annotation = columnAnnotation(df=col_anno))
dev.off()







