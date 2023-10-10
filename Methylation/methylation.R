# key steps for methylation data processing

# load package

library("ChAMP")
library("ggplot2")

# load source data file and preprocess

data <- champ.load(source.dir,arraytype='EPIC')

# filter unwanted sites

data.filter <- champ.filter()

# basic quality control

champ.QC()

# normalization with BMIQ
data.norm <- champ.norm(beta = data.filter$beta, arraytype='EPIC', cores=256)

# remove batch effect

data.remove.batch <- champ.runCombat(beta = data.norm, pd=data.filter$pd, batchname="patient.id")

# run limma through champ wrapper function

data.dmp <- champ.DMP(beta = data.norm, pheno = data.filter$pd$group)

# save out dmp

saveRDS(data.dmp,"DMP.rds")

# for starbrust plot code
# suppose you have one dataframe recording DMP and DEG

meth_sub$logpvalue.methylation <- ifelse(meth_sub$logFC>0,-log10(meth_sub$pvalue),log10(meth_sub$pvalue))

meth_sub$logpvalue.rnaseq <- ifelse(meth_sub$rnaseq.logFC>0,-log10(meth_sub$rnaseq.pvalue),log10(meth_sub$rnaseq.pvalue))

g <- ggpubr::ggscatter(meth_sub,x="logpvalue.methylation",y="logpvalue.rnaseq",color = "type")

g <- g + geom_vline(xintercept = c(log10(pvalue_cutoff),-1 * log10(pvalue_cutoff)))

g <- g + geom_hline(yintercept = c(log10(pvalue_cutoff),-1 * log10(pvalue_cutoff)))

g <- g + theme(legend.position = "bottom")

g <- g + labs(title = "Starburst Plot",x = "DNA Methylation\nLog10(Pvalue)",y = "Gene Expression\nLog10(Pvalue)")
