rm(list=ls())

# setwd("../")
library(deconstructSigs)

### read snvs
snvs=read.csv("data/0413.snv.csv",stringsAsFactors = F)

groups=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F)

pos=sapply(groups$Tumor,function(x)which(snvs$sample_id==x))

dim(snvs)
snvs=snvs[unlist(pos),]
dim(snvs)



result_dir="08_mutation_signatures"

if(!file.exists(result_dir)){
  dir.create(result_dir)
}


sigs.input <- mut.to.sigs.input(mut.ref = snvs, 
                                sample.id = "sample_id", 
                                chr = "chr", 
                                pos = "chr_start", 
                                ref = "reference_allele", 
                                alt = "variant_allele")

groups=read.csv("data/idcp.samples.with.HXid.update.csv",stringsAsFactors = F)


samples_used=unique(groups$Tumor)

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

result$sample=rownames(result)


library(reshape)


df=melt(result,id.vars = "sample")

colnames(df)=c("sample","Signature","Weight")

pos=match(df$sample,groups$Tumor)

df=cbind(df,groups[pos,])
df$group[df$group=="adenocarcinoma"]="PCA"


write.csv(df,"08_mutation_signatures/all.44samples.cosmic.signatures.csv",quote = F,row.names = F)


library(ggpubr)

df$Group=paste(df$hxid,df$group,sep="_")



p=ggbarplot(df, x="Group", y="Weight",color = "Signature",fill="Signature")+rotate_x_text(angle = 90)

print(p)
ggsave("08_mutation_signatures/all.44samples.cosmic.signatures.barplot.pdf",width = 16,height = 8)



signatures=as.character(unique(df$Signature))

my_comparisons <- list( c("IDCP", "PCA"))

pvalue_df=c()

for(cosmic_signature in signatures){
  sub=df[df$Signature==cosmic_signature,]
  pvalue=t.test(sub$Weight[sub$group=="IDCP"],sub$Weight[sub$group=="PCA"],paired = T)
  
  pvalue_df=rbind(pvalue_df,c(cosmic_signature,pvalue$p.value))
  
  if(is.na(pvalue$p.value)){
    print("NA")
  }else if(pvalue$p.value<0.05){
            pp=ggpaired(sub,x="group",y="Weight", shape = "group",color = "group",id="hxid",label = "hxid",title = cosmic_signature,
                        label.select = list(top.up = 2,top.down=2),line.color = "gray",xlab="",linetype = "twodash",ylab="Weight")+
            stat_compare_means(comparisons = my_comparisons,method = "t.test",paired = T)
            print(pp)
            ggsave(file.path("08_mutation_signatures",paste("all.44samples.cosmic",cosmic_signature,"pdf",sep=".")),width = 6,height = 6)
  }
  

}


pvalue_df=as.data.frame(pvalue_df,stringsAsFactors = F)
colnames(pvalue_df)=c("Signature","pvalue")


pvalue_df$pvalue[pvalue_df$pvalue=="NaN"]="1"
pvalue_df$pvalue=as.numeric(pvalue_df$pvalue)




ggbarplot(pvalue_df,x="Signature",y="pvalue",sort.val = "asc",rotate=TRUE,palette = "jco",
          fill="skyblue",color = "skyblue")+ylab("T test p value between IDCP and PCA")+xlab("Cosmic signatures")

ggsave("08_mutation_signatures/all.44samples.cosmic.signatures.pvalue.between.idcp.and.pca.pdf",width = 6,height = 6)
