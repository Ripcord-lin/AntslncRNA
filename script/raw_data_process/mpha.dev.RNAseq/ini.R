# cur_dir = dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(cur_dir)
library(tximport)
library(DESeq2)
library(dplyr)
tp = read.csv("Mpha_deve_sampleinfo.csv")
rownames(tp) = tp$Run
tp$Age = as.factor(tp$Age) 
tp$caste.x = as.factor(tp$caste.x)
tp$path = paste("mpha_salmon_results",tp$Run,"quant.sf",sep = "/")
info.mpha = tp
rm(tp)
# make count matrix
tx2 = read.table("tx2gene")
txi = tximport(info.mpha[["path"]],tx2gene = tx2,type = "salmon")
colnames(txi$abundance) <- info.mpha$Run
colnames(txi$counts) <- info.mpha$Run
colnames(txi$length) <- info.mpha$Run
# DESeq
Mpha_dev_dds = DESeqDataSetFromTximport(txi,colData = info.mpha,design = ~caste.x)
data.mpha = assay(vst(Mpha_dev_dds))
abundance <- txi$abundance
counts <- txi$counts
colnames(abundance) <- info.mpha$Run

geneinfo = read.csv(header = T,"Mpha_gene_info_LncAndProteinsFromChromosomes.csv")
used.geneinfo = geneinfo[which(geneinfo$gene %in% rownames(data.mpha)),]
rm(geneinfo,tx2)

library(preprocessCore)
abundance.log = log2(abundance + 1)
abundance.quantile = normalize.quantiles(as.matrix(abundance.log))
rownames(abundance.quantile) = rownames(abundance)
colnames(abundance.quantile) = colnames(abundance)
abundance.quantile <- as.data.frame(abundance.quantile)
abundance.quantile[which(abundance.quantile < 0)] = 0

print("data.mpha: vst normalized")
print("abundance: tpm")
print("counts: raw counts")
print("abundance.quantile:log2, normalize.quantiles")